#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { IQTREE } from './main.nf' 

workflow test_iqtree {

    inputs = tuple(
        [ id:"test" ],
        file(params.test_data['species']['portiera']['genome']['aln_gz'], checkIfExists: true)
    )

    IQTREE ( inputs )
}
