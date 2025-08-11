#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RGI } from './main.nf' 

workflow test_rgi {

    inputs = [
        tuple([ id:"GCF_900478275" ],
              file(params.test_data['species']['haemophilus_influenzae']['genome']['fna_gz'], checkIfExists: true)),
        tuple([ id:"GCA_000027305" ],
              file(params.test_data['species']['haemophilus_influenzae']['nthi']['fna_gz'], checkIfExists: true))
    ]

    RGI ( Channel.fromList(inputs) )
}
