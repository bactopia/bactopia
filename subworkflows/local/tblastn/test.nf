#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TBLASTN } from './main.nf' addParams(tblastn_query: params.test_data['species']['portiera']['genome']['single_protein'])
include { TBLASTN as TBLASTN_GZ } from './main.nf' addParams( tblastn_query: params.test_data['species']['portiera']['genome']['single_protein_gz'])

workflow test_tblastn {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    TBLASTN ( inputs )
}

workflow test_tblastn_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    TBLASTN_GZ ( inputs )
}
