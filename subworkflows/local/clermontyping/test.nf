#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLERMONTYPING } from './main.nf' 

workflow test_clermontyping {

    inputs = tuple(
        [ id:"GCF_001695515" ],
        file(params.test_data['species']['escherichia_coli']['genome']['fna_gz'], checkIfExists: true)
    )

    CLERMONTYPING(inputs)
}
