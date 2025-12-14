#!/usr/bin/env nextflow
/**
 * MinMER-assisted species-specific tool selection and execution.
 *
 * This Bactopia Tool, Merlin, uses MinMER distances based on the RefSeq sketch to automatically
 * run species-specific analysis tools. Merlin identifies the closest reference genomes
 * and executes appropriate typing and analysis tools for each detected species.
 *
 * @status stable
 * @keywords species-specific, automated, mash, minmer, typing, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,conditional-logic,automation
 *
 * @subworkflows bactopiatool_init, merlin
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input minmer_db
 * Path to Minmer database for species identification
 *
 * @input ask_merlin
 * Interactive mode for species-specific tool selection
 *
 * @input full_merlin
 * Execute all species-specific tools regardless of species match
 *
 * @section Species-Specific Analysis
 * @note Tools executed depend on detected species
 * @publish                         Analysis results from all executed species-specific tools
 *
 * @section Merged Results
 * @publish merlin.tsv                Merged summary of all species-specific analyses
 *
 * @section Execution Logs
 * @publish logs/**                   Tool execution logs from all executed tools
 * @publish logs/nf-*                 Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml              Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    emmtyper_blastdb      : Path?
    hicap_database_dir    : Path?
    hicap_model_fp        : Path?
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
include { MERLIN            } from '../../../subworkflows/merlin/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    DATASETS()
    MERLIN(
        BACTOPIATOOL_INIT.out.assembly_reads,
        DATASETS.out.mash_db,
        // emmtyper
        params.emmtyper_blastdb,
        // hicap
        params.hicap_database_dir,
        params.hicap_model_fp,
        // staphtyper
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )

    // Collect outputs
    ch_results = ch_results.mix(MERLIN.out.results)
    ch_logs = ch_logs.mix(MERLIN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MERLIN.out.nf_logs)
    ch_versions = ch_versions.mix(MERLIN.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
