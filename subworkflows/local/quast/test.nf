#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QUAST } from './main.nf' 

workflow test_quast {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['meta'], checkIfExists: true)
    )

    QUAST ( inputs )
}
