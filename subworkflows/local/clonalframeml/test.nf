#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLONALFRAMEML } from './main.nf' 

workflow test_clonalframeml {

    input = tuple(
        [ id:'test' ], // meta map
        file(params.test_data['species']['portiera']['genome']['aln_gz'], checkIfExists: true)
    )

    CLONALFRAMEML(input)
}
