#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EMMTYPER } from './main.nf' 

workflow test_emmtyper {

    inputs = tuple( )

    EMMTYPER ( inputs )
}
