#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { ABRICATE          } from '../../../subworkflows/abricate/main'
include { workflowSummary   } from 'plugin/nf-bactopia'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
//WorkflowMain.initialise(workflow, params, log)
//WorkflowBactopiaTools.initialise(workflow, params, log)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    BACTOPIATOOL_INIT(params.validate_params)
    ABRICATE(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = ABRICATE.out.tsv.mix(ABRICATE.out.merged_tsv)
    logs = ABRICATE.out.logs
    nf_logs = ABRICATE.out.nf_logs
    versions = ABRICATE.out.versions
}



/*
    tsv = ABRICATE_RUN.out.report
    run_logs = ABRICATE_RUN.out.logs
    run_nf_logs = ABRICATE_RUN.out.nf_logs
    summary_logs = ABRICATE_SUMMARY.out.logs
    summary_nf_logs = ABRICATE_SUMMARY.out.nf_logs
    merged_tsv = ch_merged_abricate
    versions = ch_versions // channel: [ versions.yml ]
*/

output {
    'results' {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    'logs' {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    'nf_logs' {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    'versions' {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}


/*
========================================================================================
    THE END
========================================================================================
*/
