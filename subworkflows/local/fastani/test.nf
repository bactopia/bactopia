#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTANI } from './main.nf' 

workflow test_fastani {

    inputs = [
        tuple([ id:"GCF_900478275" ],
              file(params.test_data['species']['haemophilus_influenzae']['genome']['fna_gz'], checkIfExists: true)),
        tuple([ id:"GCA_000027305" ],
              file(params.test_data['species']['haemophilus_influenzae']['nthi']['fna_gz'], checkIfExists: true))
    ]

    reference = [
        file(params.test_data['species']['haemophilus_influenzae']['genome']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['haemophilus_influenzae']['nthi']['fna_gz'], checkIfExists: true)
    ]

    FASTANI(Channel.fromList(inputs), Channel.fromList(reference))
}
