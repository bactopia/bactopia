#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SSUISSERO } from './main.nf' 

workflow test_ssuissero {

    inputs = tuple(
        [ id:"GCF_002285535" ],
        file(params.test_data['species']['streptococcus_suis']['genome']['fna_gz'], checkIfExists: true)
    )

    SSUISSERO(inputs)
}
