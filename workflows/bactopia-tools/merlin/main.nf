#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
include { MERLIN            } from '../../../subworkflows/merlin/main'
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
    DATASETS()
    MERLIN(
        BACTOPIATOOL_INIT.out.samples,
        DATASETS.out.mash_db,
        // emmtyper
        params.emmtyper_blastdb ? file(params.emmtyper_blastdb, checkIfExists: true) : [],
        // hicap
        params.hicap_database_dir ? file(params.hicap_database_dir, checkIfExists: true) : [],
        params.hicap_model_fp ? file(params.hicap_model_fp, checkIfExists: true) : [],
        // staphtyper
        params.spatyper_repeats ? file(params.spatyper_repeats, checkIfExists: true) : [],
        params.spatyper_repeat_order ? file(params.spatyper_repeat_order, checkIfExists: true) : []
    )

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = MERLIN.out.results
    logs = MERLIN.out.logs
    nf_logs = MERLIN.out.nf_logs
    versions = MERLIN.out.versions
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
