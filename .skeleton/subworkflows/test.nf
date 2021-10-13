#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SUBWORKFLOW } from './main.nf' 

workflow test_SUBWORKFLOW {

    inputs = tuple( )

    SUBWORKFLOW ( inputs )
}
