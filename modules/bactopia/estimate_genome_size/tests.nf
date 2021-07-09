#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ESTIMATE_GENOME_SIZE } from './main.nf' 

workflow test_estimate_genome_size_pe {

    inputs = tuple(
        "test_estimate_genome_size_pe",
        "paired-end",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ESTIMATE_GENOME_SIZE ( inputs, 0 )
}

workflow test_estimate_genome_size_pe_known {

    inputs = tuple(
        "test_estimate_genome_size_pe_known",
        "paired-end",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ESTIMATE_GENOME_SIZE ( inputs, params.test_data['reference']['length'] )
}


workflow test_estimate_genome_size_se {

    inputs = tuple(
        "test_estimate_genome_size_se",
        "single-end",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ESTIMATE_GENOME_SIZE ( inputs, 0 )
}
