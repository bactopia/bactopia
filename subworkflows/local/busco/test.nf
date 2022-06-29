#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BUSCO } from './main.nf' 

workflow test_busco {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    BUSCO ( inputs )
}
