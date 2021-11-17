#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BAKTA } from './main.nf' 

workflow test_bakta {

    inputs = tuple( )

    BAKTA ( inputs )
}
