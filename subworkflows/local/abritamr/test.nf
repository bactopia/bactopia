#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ABRITAMR } from './main.nf' 

workflow test_abritamr {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['escherichia_coli_mcr1']['genome']['fna_gz'], checkIfExists: true)
    )

    ABRITAMR ( inputs )
}
