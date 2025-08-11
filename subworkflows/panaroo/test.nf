#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PANAROO } from './main.nf' 

workflow test_panaroo {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['species']['portiera']['genome']['gff_gz'], checkIfExists: true),
                file(params.test_data['species']['portiera']['genome']['gff2_gz'], checkIfExists: true),
                file(params.test_data['species']['portiera']['genome']['gff3_gz'], checkIfExists: true) ]
    ]

    PANAROO(input)
}
