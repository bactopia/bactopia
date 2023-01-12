#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TETON } from './main.nf' 

workflow test_teton {

    inputs = tuple( )

    TETON ( inputs )
}
