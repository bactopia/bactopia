#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_READS } from './main.nf' 

workflow test_qc_reads_pe {

    inputs = tuple(
        "test_qc_reads_pe",
        "paired-end",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['reference']['genome_size'], checkIfExists: true)
    )

    QC_READS ( inputs )
}

workflow test_qc_reads_se {

    inputs = tuple(
        "test_qc_reads_se",
        "single-end",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['reference']['genome_size'], checkIfExists: true)
    )

    QC_READS ( inputs )
}
