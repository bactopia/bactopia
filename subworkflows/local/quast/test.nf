#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QUAST } from './main.nf' 

workflow test_quast {

    inputs = tuple( )

    QUAST ( inputs )
}
