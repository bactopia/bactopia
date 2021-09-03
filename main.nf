#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

//include { BACTOPIA } from './workflows/bactopia'
//include { STAPHOPIA } from './workflows/staphopia'

//
// WORKFLOW: Run main nf-core/bactmap analysis pipeline
//
//workflow NFCORE_BACTMAP {
//    BACTOPIA ()
//}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
//workflow {
//    NFCORE_BACTMAP ()
//}

/*
========================================================================================
    THE END
========================================================================================
*/
