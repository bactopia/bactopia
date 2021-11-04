#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NGMASTER } from './main.nf' 

workflow test_ngmaster {

    inputs = tuple(
        [ id:"GCF_001047255" ],
        file(params.test_data['species']['neisseria_gonorrhoeae']['genome']['fna_gz'], checkIfExists: true)
    )

    NGMASTER ( inputs )
}
