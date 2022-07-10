#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MCRONI } from './main.nf' 

workflow test_mcroni {

    inputs = tuple(
        [ id:"GCF_001682305" ],
        file(params.test_data['species']['escherichia_coli_mcr1']['genome']['fna_gz'], checkIfExists: true)
    )

    MCRONI ( inputs )
}
