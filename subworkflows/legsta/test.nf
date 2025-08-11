#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { LEGSTA } from './main.nf' 

workflow test_legsta {

    inputs = tuple(
        [ id:"GCF_000048645" ],
        file(params.test_data['species']['legionella_pneumophila']['genome']['fna_gz'], checkIfExists: true)
    )

    LEGSTA ( inputs )
}
