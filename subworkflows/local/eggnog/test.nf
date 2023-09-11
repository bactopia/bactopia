#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EGGNOG } from './main.nf' addParams( eggnog_db: params.test_data['datasets']['eggnog'] )
include { EGGNOG as EGGNOG_TARBALL } from './main.nf' addParams( eggnog_db: params.test_data['datasets']['eggnog_tarball'] )

workflow test_eggnog {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['faa_single_gz'], checkIfExists: true)
    )

    EGGNOG ( inputs )
}

workflow test_eggnog_tarball {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['faa_single_gz'], checkIfExists: true)
    )

    EGGNOG_TARBALL ( inputs )
}

