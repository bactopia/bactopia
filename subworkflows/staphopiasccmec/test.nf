#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STAPHOPIASCCMEC } from './main.nf' 

workflow test_staphopiasccmec {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna'], checkIfExists: true)
    )

    STAPHOPIASCCMEC ( inputs )
}
