#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SHIGAPASS } from './main.nf' 

workflow test_shigapass {

    inputs = tuple(
        [ id:"GCF_016726285" ],
        file(params.test_data['species']['shigella_boydii']['genome']['fna_gz'], checkIfExists: true)
    )

    SHIGAPASS ( inputs )
}
