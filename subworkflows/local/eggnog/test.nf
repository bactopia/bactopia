#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EGGNOG } from './main.nf' addParams( eggnog: params.test_data['datasets']['eggnog'], download_eggnog : true )

workflow test_eggnog {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['faa_single_gz'], checkIfExists: true)
    )

    EGGNOG ( inputs )
}
