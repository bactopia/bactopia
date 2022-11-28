#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PASTY } from './main.nf' 

workflow test_pasty {

    inputs = tuple(
        [ id:"GCF_000006765" ],
        file(params.test_data['species']['pseudomonas_aeruginosa']['genome']['fna_gz'], checkIfExists: true)
    )

    PASTY ( inputs )
}


