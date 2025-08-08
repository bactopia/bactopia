#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ISMAPPER_WORKFLOW } from './main.nf'

workflow test_ismapper {
    ch_reads = Channel.fromPath(params.test_data['species']['portiera']['illumina']['r{1,2}'], checkIfExists: true)
        .map { files -> [[id: 'test'], files] }
    
    ISMAPPER_WORKFLOW(
        ch_reads,
        params.test_data['species']['portiera']['genome']['gbk'],
        params.test_data['species']['portiera']['insertion']['fasta']
    )
}
