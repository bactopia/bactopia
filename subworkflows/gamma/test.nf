#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GAMMA } from './main.nf' addParams(gamma_db: params.test_data['species']['portiera']['genome']['ffn'])

workflow test_gamma {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    GAMMA ( inputs )
}
