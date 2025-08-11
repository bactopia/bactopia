#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLASTP } from './main.nf' addParams(blastp_query: params.test_data['species']['portiera']['genome']['single_protein'])
include { BLASTP as BLASTP_GZ } from './main.nf' addParams( blastp_query: params.test_data['species']['portiera']['genome']['single_protein_gz'])
include { BLASTP as BLASTP_USE_GENES } from './main.nf' addParams( blastp_use_genes: true, blastp_query: params.test_data['species']['portiera']['genome']['single_protein'])

workflow test_blastp {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTP ( inputs )
}

workflow test_blastp_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTP_GZ ( inputs )
}

workflow test_blastp_use_genes {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTP_USE_GENES ( inputs )
}
