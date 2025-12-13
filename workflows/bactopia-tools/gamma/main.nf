#!/usr/bin/env nextflow
/**
 * Identification, classification, and annotation of translated gene matches.
 *
 * This Bactopia Tool uses [GAMMA](https://github.com/rastanton/GAMMA) to identify, classify, and annotate
 * translated gene matches from assemblies using a comprehensive protein database.
 *
 * @status stable
 * @keywords gene annotation, protein classification, translation, gamma, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,annotation
 * @citation gamma
 *
 * @subworkflows bactopiatool_init, gamma
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input gamma_db
 * Path to GAMMA database
 *
 * @section Per-Sample Results
 * @publish *.tsv            Tab-delimited file with gene classification results
 *
 * @section Merged Results
 * @publish gamma.tsv         Merged TSV file containing GAMMA results from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    gamma_db : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { GAMMA             } from '../../../subworkflows/gamma/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    GAMMA(
        BACTOPIATOOL_INIT.out.samples,
        params.gamma_db
    )

    // Collect outputs
    ch_results = ch_results.mix(GAMMA.out.results)
    ch_logs = ch_logs.mix(GAMMA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GAMMA.out.nf_logs)
    ch_versions = ch_versions.mix(GAMMA.out.versions)

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
