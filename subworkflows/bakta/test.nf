#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.bakta_db = "${params.test_data['datasets']['bakta']}"
include { BAKTA } from './main.nf'
include { BAKTA as BAKTA_TARBALL } from './main.nf'

workflow test_bakta {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    BAKTA ( inputs )
}

workflow test_bakta_tarball {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    BAKTA_TARBALL ( inputs )
}
