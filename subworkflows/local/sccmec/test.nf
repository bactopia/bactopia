#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCCMEC } from './main.nf' 

workflow test_sccmec {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna'], checkIfExists: true)
    )

    SCCMEC ( inputs )
}


