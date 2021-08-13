#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ASSEMBLE_GENOME } from './main.nf' 

workflow test_assemble_genome_se {
    // sample, sample_type, single_end, fastqs, extra, genome_size
    inputs = tuple(
        "test_assemble_genome_se",
        'single-end',
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna']),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_pe {
    // sample, sample_type, single_end, fastqs, extra, genome_size
    inputs = tuple(
        "test_assemble_genome_pe",
        'paired-end',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna']),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_hybrid {
    // sample, sample_type, single_end, fastqs, extra, genome_size
    inputs = tuple(
        "test_assemble_genome_hybrid",
        'hybrid',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['ont']['se'], checkIfExists: true),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}