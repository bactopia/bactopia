#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PBPTYPER } from './main.nf' 

workflow test_pbptyper {

    inputs = tuple(
        [ id:"GCF_001457635" ],
        file(params.test_data['species']['streptococcus_pneumoniae']['genome']['fna_gz'], checkIfExists: true)
    )

    PBPTYPER ( inputs )
}


