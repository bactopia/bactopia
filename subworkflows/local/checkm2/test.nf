#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CHECKM2 } from './main.nf' 

workflow test_checkm2 {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    CHECKM2 ( inputs )
}
