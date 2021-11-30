#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KLEBORATE } from './main.nf' 

workflow test_kleborate {

    inputs = tuple(
        [ id:"GCF_000009885" ],
        file(params.test_data['species']['klebsiella_pneumoniae']['genome']['fna_gz'], checkIfExists: true)
    )

    KLEBORATE(inputs)
}
