#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ABRITAMR } from './main.nf' 

workflow test_abritamr {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['faa_gz'], checkIfExists: true)
    )

    ABRITAMR ( inputs )
}
