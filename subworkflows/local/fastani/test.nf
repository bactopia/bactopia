#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTANI } from './main.nf' 

workflow test_fastani {

    inputs = tuple( )

    FASTANI ( inputs )
}
