#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLASTX } from './main.nf' addParams(blastx_query: params.test_data['species']['portiera']['genome']['single_gene'])
include { BLASTX as BLASTX_GZ } from './main.nf' addParams( blastx_query: params.test_data['species']['portiera']['genome']['single_gene_gz'])
include { BLASTX as BLASTX_USE_GENES } from './main.nf' addParams( blastx_use_genes: true, blastx_query: params.test_data['species']['portiera']['genome']['single_gene'])

workflow test_blastx {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTX ( inputs )
}

workflow test_blastx_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTX_GZ ( inputs )
}

workflow test_blastx_use_genes {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTX_USE_GENES ( inputs )
}
