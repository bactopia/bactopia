#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PANGENOME } from './main.nf'

workflow test_pangenome {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['species']['portiera']['genome']['gff_gz'], checkIfExists: true),
                file(params.test_data['species']['portiera']['genome']['gff2_gz'], checkIfExists: true),
                file(params.test_data['species']['portiera']['genome']['gff3_gz'], checkIfExists: true),
                file(params.test_data['species']['portiera']['genome']['gff4_gz'], checkIfExists: true) ]
    ]

    PANGENOME ( input )
}
