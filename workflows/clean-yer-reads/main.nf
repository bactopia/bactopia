#!/usr/bin/env nextflow
nextflow.preview.output = true

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
    // Check if help is requested
    if (params.help || params.help_all) {
        log.info paramsHelp()
        exit 0
    }

    // Initialize and execute the workflow
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()
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

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = ch_results
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
