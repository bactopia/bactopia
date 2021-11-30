#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { LISSERO } from './main.nf' 

workflow test_lissero {

    inputs = tuple(
        [ id:"GCF_002285835" ],
        file(params.test_data['species']['listeria_monocytogenes']['genome']['fna_gz'], checkIfExists: true)
    )

    LISSERO ( inputs )
}
