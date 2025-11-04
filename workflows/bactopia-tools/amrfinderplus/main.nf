#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { AMRFINDERPLUS     } from '../../../subworkflows/amrfinderplus/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
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

    if (params.amrfinder_db) {
        // User specified database
        AMRFINDERPLUS(BACTOPIATOOL_INIT.out.samples, file(params.amrfinder_db))
    } else {
        // Use default database
        DATASETS()
        AMRFINDERPLUS(BACTOPIATOOL_INIT.out.samples, DATASETS.out.amrfinderplus_db)
    }

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = AMRFINDERPLUS.out.results
    logs = AMRFINDERPLUS.out.logs
    nf_logs = AMRFINDERPLUS.out.nf_logs
    versions = AMRFINDERPLUS.out.versions
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
