#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATHER_FASTQS } from './main.nf' 

workflow test_gather_fastqs_pe {

    inputs = tuple(
        'test_gather_fastqs_pe',
        'paired-end',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_merge_pe {

    inputs = tuple(
        'test_gather_fastqs_merge_pe',
        'merge-pe',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_se {

    inputs = tuple(
        'test_gather_fastqs_se',
        'single-end',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_merge_se {

    inputs = tuple(
        'test_gather_fastqs_merge_se',
        'merge-se',
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_sra_accession {

    inputs = tuple(
        params.test_data['accessions']['srx'],
        'sra-accession',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_assembly_accession {

    inputs = tuple(
        params.test_data['accessions']['gcf'],
        'assembly-accession',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}

workflow test_gather_fastqs_assembly {

    inputs = tuple(
        params.test_data['reference']['name'],
        'assembly',
        false,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    GATHER_FASTQS ( inputs )
}