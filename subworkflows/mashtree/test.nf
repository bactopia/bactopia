#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASHTREE } from './main.nf' 

workflow test_mashtree {

    inputs = tuple(
        [ id:"mashtree" ],
        [file(params.test_data['species']['klebsiella_pneumoniae']['genome']['fna_gz'], checkIfExists: true),
         file(params.test_data['species']['staphylococcus_aureus']['genome']['fna_gz'], checkIfExists: true)]
    )

    MASHTREE ( inputs )
}
