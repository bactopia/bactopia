#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PROKKA } from './main.nf' 

workflow test_prokka {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    PROKKA ( inputs )
}
