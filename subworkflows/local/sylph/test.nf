#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SYLPH } from './main.nf' 

workflow test_sylph {

    inputs = tuple( )

    SYLPH ( inputs )
}
