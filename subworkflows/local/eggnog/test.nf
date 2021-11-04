#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EGGNOG } from './main.nf' 

workflow test_eggnog {

    inputs = tuple( )

    EGGNOG ( inputs )
}
