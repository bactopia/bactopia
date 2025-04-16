#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'

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
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    BACTOPIATOOL_INIT (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

}

/*
========================================================================================
    THE END
========================================================================================
*/
