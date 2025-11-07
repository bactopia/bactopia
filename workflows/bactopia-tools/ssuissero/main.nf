#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { SSUISSERO     } from '../../../subworkflows/ssuissero/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Initialize and execute the workflow
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()

    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    SSUISSERO(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(SSUISSERO.out.results)
    ch_logs = ch_logs.mix(SSUISSERO.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SSUISSERO.out.nf_logs)
    ch_versions = ch_versions.mix(SSUISSERO.out.versions)

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
    run_results {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs {
        path { meta, file -> {
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        } }
    }
    run_versions {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    sample_versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
