#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DEFENSEFINDER } from './main.nf' 

workflow test_defensefinder {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    DEFENSEFINDER ( inputs )
}
