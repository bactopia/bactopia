#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SISTR } from './main.nf' 

workflow test_sistr {

    inputs = tuple(
        [ id:"GCF_016028495" ],
        file(params.test_data['species']['salmonella_enterica']['genome']['fna_gz'], checkIfExists: true)
    )

    SISTR ( inputs )
}
