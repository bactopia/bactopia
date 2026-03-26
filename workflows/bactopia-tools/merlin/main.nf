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
include { DATASETS          } from '../../../subworkflows/bactopia/datasets/main'
include { MERLIN            } from '../../../subworkflows/merlin/main'

workflow {
    main:
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

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = MERLIN.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MERLIN.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = MERLIN.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = MERLIN.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
