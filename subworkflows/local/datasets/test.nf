#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DATASETS } from './main.nf' 

workflow test_datasets {

    inputs = tuple( )

    DATASETS ( inputs )
}
