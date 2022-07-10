#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GUBBINS } from './main.nf' 

workflow test_gubbins {

    input = tuple(
        [ id:'test' ], // meta map
        file(params.test_data['species']['portiera']['genome']['aln_gz'], checkIfExists: true)
    )

    GUBBINS(input)
}
