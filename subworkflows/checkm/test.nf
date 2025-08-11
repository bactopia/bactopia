#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CHECKM } from './main.nf' 

workflow test_checkm {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    CHECKM ( inputs )
}
