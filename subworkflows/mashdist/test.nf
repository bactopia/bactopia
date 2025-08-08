#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASHDIST } from './main.nf'

workflow test_mashdist {
    ch_seqs = Channel.fromPath(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
        .map { genome -> [[id: genome.baseName], genome] }
    
    MASHDIST(ch_seqs, params.test_data['species']['portiera']['sketches']['sketch'])
}
