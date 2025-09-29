#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { STAPHTYPER        } from '../../../subworkflows/staphtyper/main'
include { paramsHelp        } from 'plugin/nf-bactopia'
include { workflowSummary   } from 'plugin/nf-bactopia'

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
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    STAPHTYPER(
        BACTOPIATOOL_INIT.out.samples,
        params.repeats ? file(params.repeats, checkIfExists: true) : [],
        params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []
    )

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = STAPHTYPER.out.results
    logs = STAPHTYPER.out.logs
    nf_logs = STAPHTYPER.out.nf_logs
    versions = STAPHTYPER.out.versions
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
