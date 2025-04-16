#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ABRICATE } from './main.nf' 

workflow test_abricate {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    ABRICATE ( inputs )
}
