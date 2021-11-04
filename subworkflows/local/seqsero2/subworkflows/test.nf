#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SUBWORKFLOW_UPPER } from './main.nf' 

workflow test_SUBWORKFLOW_NAME {

    inputs = tuple( )

    SUBWORKFLOW_UPPER ( inputs )
}
