#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STAPHTYPER } from './main.nf' 

workflow test_staphtyper {

    inputs = tuple(
        [ id:"output" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna'], checkIfExists: true)
    )

    STAPHTYPER ( inputs )
}
