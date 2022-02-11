#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { AMRFINDERPLUS } from './main.nf' 

workflow test_amrfinderplus {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['faa_gz'], checkIfExists: true)
    )

    AMRFINDERPLUS ( inputs )
}
