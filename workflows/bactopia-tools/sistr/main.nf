#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { SISTR         } from '../../../subworkflows/sistr/main'
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
    SISTR(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = SISTR.out.results
    logs = SISTR.out.logs
    nf_logs = SISTR.out.nf_logs
    versions = SISTR.out.versions
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
