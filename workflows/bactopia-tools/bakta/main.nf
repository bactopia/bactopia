#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { BAKTA             } from '../../../subworkflows/bakta/main'
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
    BAKTA(
        BACTOPIATOOL_INIT.out.samples,
        params.bakta_db ? file(params.bakta_db) : [],
        params.download_bakta,
        params.bakta_save_as_tarball,
        params.bakta_proteins ? file(params.bakta_proteins, checkIfExists: true) : [],
        params.bakta_prodigal_tf ? file(params.bakta_prodigal_tf, checkIfExists: true) : [],
        params.replicons ? file(params.replicons, checkIfExists: true) : []
    )

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = BAKTA.out.results
    logs = BAKTA.out.logs
    nf_logs = BAKTA.out.nf_logs
    versions = BAKTA.out.versions
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
