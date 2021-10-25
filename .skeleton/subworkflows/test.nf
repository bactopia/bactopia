#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SUBWORKFLOW_NAME } from './main.nf' 

workflow test_SUBWORKFLOW_NAME {

    inputs = tuple( )

    SUBWORKFLOW_NAME ( inputs )
}
