#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQSERO2 } from './main.nf'

workflow test_seqsero2 {

    inputs = tuple(
        [ id:"GCF_016028495" ],
        file(params.test_data['species']['salmonella_enterica']['genome']['fna_gz'], checkIfExists: true)
    )

    SEQSERO2 ( inputs )
}
