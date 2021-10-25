#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (params.workflows.containsKey(params.wf)) {
    if (params.workflows[params.wf].containsKey('is_workflow')) {
        if (params.wf == "staphopia") {
            include { STAPHOPIA } from './workflows/staphopia'
        } else {
            include { BACTOPIA } from './workflows/bactopia'
        }
    } else {
        include { BACTOPIATOOLS } from './workflows/bactopia-tools'
    }
} else {
    log.error "${params.wf} is not an available Bactopia Tool. Use --available_tools to see full list"
    exit 1
}

/*
========================================================================================
    RUN WORKFLOWS
========================================================================================
*/

workflow {
    if (params.workflows[params.wf].containsKey('is_workflow')) {
        if (params.wf == "staphopia") {
            STAPHOPIA()
        } else {
            BACTOPIA()
        }
    } else {
        BACTOPIATOOLS()
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
