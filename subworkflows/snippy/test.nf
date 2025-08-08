#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SNIPPY } from './main.nf'

workflow test_snippy {
    ch_reads = Channel.fromPath(params.test_data['species']['portiera']['illumina']['r{1,2}'], checkIfExists: true)
        .map { files -> [[id: 'test', single_end: false, runtype: 'paired-end'], files] }
    
    ch_reference = Channel.fromPath(params.test_data['species']['portiera']['genome']['gbk'], checkIfExists: true)
        .map { reference -> [[id: 'reference'], reference] }
    
    SNIPPY(ch_reads, ch_reference)
}
