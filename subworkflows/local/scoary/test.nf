#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCOARY } from './main.nf' 

workflow test_scoary {

    inputs = tuple( )

    SCOARY ( inputs )
}
