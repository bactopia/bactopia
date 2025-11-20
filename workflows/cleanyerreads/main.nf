#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Core
include { BACTOPIA_INIT   } from '../../subworkflows/utils/bactopia'
include { GATHER          } from '../../subworkflows/bactopia/gather/main'
include { QC              } from '../../subworkflows/bactopia/qc/main'
include { paramsHelp      } from 'plugin/nf-bactopia'
include { workflowSummary } from 'plugin/nf-bactopia'

// Scrubber
include { SCRUBBER        } from '../../subworkflows/scrubber/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {
    main:
    // Initialize and execute the workflow
    ch_results = channel.empty()
    ch_logs = channel.empty()
    ch_nf_logs = channel.empty()
    ch_versions = channel.empty()
    BACTOPIA_INIT()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)
    ch_results = ch_results.mix(GATHER.out.results)
    ch_logs = ch_logs.mix(GATHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
    ch_versions = ch_versions.mix(GATHER.out.versions)

    if (params.use_k2scrubber || params.use_srascrubber) {
        // Remove host reads
        SCRUBBER(GATHER.out.fastq_only, params.use_srascrubber)
        ch_results = ch_results.mix(SCRUBBER.out.results)
        ch_logs = ch_logs.mix(SCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SCRUBBER.out.nf_logs)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Clean up scrubbed reads
        QC(
            SCRUBBER.out.scrubbed_extra,
            params.adapters ? file(params.adapters, checkIfExists: true) : [],
            params.phix ? file(params.phix, checkIfExists: true) : []
        )
    } else {
        // Clean up raw reads
        QC(
            GATHER.out.raw_fastq,
            params.adapters ? file(params.adapters, checkIfExists: true) : [],
            params.phix ? file(params.phix, checkIfExists: true) : []
        )
    }
    ch_results = ch_results.mix(QC.out.results)
    ch_logs = ch_logs.mix(QC.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QC.out.nf_logs)
    ch_versions = ch_versions.mix(QC.out.versions)

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
