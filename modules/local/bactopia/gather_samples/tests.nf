#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATHER_SAMPLES } from './main.nf' 

workflow test_gather_samples_pe {

    inputs = tuple(
        'test_gather_samples_pe',
        'paired-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_pe {

    inputs = tuple(
        'test_gather_samples_merge_pe',
        'merge-pe',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_pe_hybrid {

    inputs = tuple(
        'test_gather_samples_merge_pe_hybrid',
        'hybrid-merge-pe',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['nanopore']['se'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_se {

    inputs = tuple(
        'test_gather_samples_se',
        'single-end',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_se {

    inputs = tuple(
        'test_gather_samples_merge_se',
        'merge-se',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_sra_accession {

    inputs = tuple(
        params.test_data['accessions']['srx'],
        'sra-accession',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_assembly_accession {

    inputs = tuple(
        params.test_data['accessions']['gcf'],
        'assembly-accession',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_assembly {

    inputs = tuple(
        params.test_data['reference']['name'],
        'assembly',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}
