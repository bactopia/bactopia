#!/usr/bin/env nextflow
/**
 * Bactopia Tool: Amrfinderplus.
 *
 * Identify antimicrobial resistance genes and point mutations in bacterial genomes.
 * This Bactopia Tool uses [AMRFinder+](https://github.com/ncbi/amr) to screen assemblies and proteins
 * for antimicrobial resistance genes, virulence genes, and resistance-associated point mutations.
 * It identifies acquired AMR genes and some point mutations in protein or assembled nucleotide sequences.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance, virulence, genes, proteins, mutations
 * @tags complexity:moderate input-type:directory output-type:multiple features:database-dependent,amr-detection
 * @citation amrfinderplus
 *
 * @subworkflows bactopiatool_init, amrfinderplus
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @input amrfinder_db
 * Path to AMRFinder+ database for AMR gene detection
 *
 * @section Per-Sample Results
 * @publish *.tsv                    AMR gene detection results in TSV format
 * @publish *-mutations.tsv          Point mutations associated with antimicrobial resistance
 *
 * @section Merged Results
 * @publish amrfinderplus.tsv        Combined AMR detection results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                  Tool execution logs
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml              Software version information
   */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    amrfinder_db : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { AMRFINDERPLUS     } from '../../../subworkflows/amrfinderplus/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()

    if (params.amrfinder_db) {
        // User specified database
        AMRFINDERPLUS(
            BACTOPIATOOL_INIT.out.assembly_proteins_gff,
            params.amrfinder_db
        )
    } else {
        // Use default database
        DATASETS()
        AMRFINDERPLUS(
            BACTOPIATOOL_INIT.out.assembly_proteins_gff,
            DATASETS.out.amrfinderplus_db
        )
    }

    // Collect outputs
    ch_results = ch_results.mix(AMRFINDERPLUS.out.results)
    ch_logs = ch_logs.mix(AMRFINDERPLUS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(AMRFINDERPLUS.out.nf_logs)
    ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions)

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
