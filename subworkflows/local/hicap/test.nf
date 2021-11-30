#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HICAP } from './main.nf' 

workflow test_hicap {

    inputs = [
        tuple([ id:"GCF_900478275" ],
              file(params.test_data['species']['haemophilus_influenzae']['genome']['fna_gz'], checkIfExists: true)),
        tuple([ id:"GCA_000027305" ],
              file(params.test_data['species']['haemophilus_influenzae']['nthi']['fna_gz'], checkIfExists: true))
    ]

    HICAP(Channel.fromList(inputs))
}
