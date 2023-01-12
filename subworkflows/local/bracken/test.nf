#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BRACKEN } from './main.nf' 

workflow test_bracken {

    inputs = tuple( )

    BRACKEN ( inputs )
}
