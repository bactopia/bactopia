#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NGMASTER } from './main.nf' 

workflow test_ngmaster {

    inputs = tuple( )

    NGMASTER ( inputs )
}
