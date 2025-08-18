#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { ISMAPPER_WORKFLOW } from '../../../subworkflows/ismapper/main'
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
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }

    // Initialize and execute the workflow
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    
    // Reference and insertion sequences should be provided via params
    ch_reference = file(params.reference, checkIfExists: true)
    ch_insertions = file(params.insertions, checkIfExists: true)
    ISMAPPER_WORKFLOW(BACTOPIATOOL_INIT.out.samples, ch_reference, ch_insertions)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = ISMAPPER_WORKFLOW.out.results
    logs = ISMAPPER_WORKFLOW.out.logs
    nf_logs = ISMAPPER_WORKFLOW.out.nf_logs
    versions = ISMAPPER_WORKFLOW.out.versions
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
