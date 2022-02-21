#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HPSUISSERO } from './main.nf' 

workflow test_hpsuissero {

    inputs = tuple(
        [ id:"GCF_002777395" ],
        file(params.test_data['species']['glaesserella_parasuis']['genome']['fna_gz'], checkIfExists: true)
    )

    HPSUISSERO(inputs)
}
