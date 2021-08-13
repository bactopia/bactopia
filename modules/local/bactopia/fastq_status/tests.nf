#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQ_STATUS } from './main.nf' 

workflow test_fastq_status_pe {

    inputs = tuple(
        "test_fastq_status_pe",
        "paired-end",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    FASTQ_STATUS ( inputs )
}

workflow test_fastq_status_se {

    inputs = tuple(
        "test_fastq_status_se",
        "single-end",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    FASTQ_STATUS ( inputs )
}
