#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATHER } from './main.nf' 

workflow test_gather {

    inputs = tuple( )

    GATHER ( inputs )
}
