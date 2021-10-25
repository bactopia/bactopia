#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SPATYPER } from './main.nf' 

workflow test_spatyper {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna'], checkIfExists: true)
    )

    SPATYPER ( inputs )
}
