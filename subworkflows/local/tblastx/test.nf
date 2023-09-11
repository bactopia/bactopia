#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TBLASTX } from './main.nf' addParams(tblastx_query: params.test_data['species']['portiera']['genome']['single_gene'])
include { TBLASTX as TBLASTX_GZ } from './main.nf' addParams( tblastx_query: params.test_data['species']['portiera']['genome']['single_gene_gz'])
include { TBLASTX as TBLASTX_USE_GENES } from './main.nf' addParams( tblastx_use_genes: true, tblastx_query: params.test_data['species']['portiera']['genome']['single_gene'])

workflow test_tblastx {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    TBLASTX ( inputs )
}

workflow test_tblastx_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    TBLASTX_GZ ( inputs )
}
