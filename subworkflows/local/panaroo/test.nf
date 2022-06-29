#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PANAROO } from './main.nf' 

workflow test_panaroo {

    inputs = tuple( )

    PANAROO ( inputs )
}
