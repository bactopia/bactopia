#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLASTN } from './main.nf' addParams(blastn_query: params.test_data['species']['portiera']['genome']['single_gene'])
include { BLASTN as BLASTN_GZ } from './main.nf' addParams( blastn_query: params.test_data['species']['portiera']['genome']['single_gene_gz'])
include { BLASTN as BLASTN_USE_GENES } from './main.nf' addParams( blastn_use_genes: true, blastn_query: params.test_data['species']['portiera']['genome']['single_gene'])

workflow test_blastn {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN ( inputs )
}

workflow test_blastn_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN_GZ ( inputs )
}

workflow test_blastn_use_genes {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN_USE_GENES ( inputs )
}
