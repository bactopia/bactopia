#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PLASMIDFINDER } from './main.nf' 

workflow test_plasmidfinder {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna_gz'], checkIfExists: true)
    )

    PLASMIDFINDER ( inputs )
}
