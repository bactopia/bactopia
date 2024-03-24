#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SRAHUMANSCRUBBER } from './main.nf' 

workflow test_srahumanscrubber {

    inputs = tuple( )

    SRAHUMANSCRUBBER ( inputs )
}
