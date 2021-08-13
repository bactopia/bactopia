#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DOWNLOAD_REFERENCES } from './main.nf' 

workflow test_download_references_pe {

    inputs = tuple(
        "test_download_references_pe",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['illumina']['msh'], checkIfExists: true)
    )

    DOWNLOAD_REFERENCES ( inputs, file(params.test_data['reference']['msh'], checkIfExists: true) )
}

workflow test_download_references_se {

    inputs = tuple(
        "test_download_references_se",
        true,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        file(params.test_data['illumina']['msh'], checkIfExists: true)
    )

    DOWNLOAD_REFERENCES ( inputs, file(params.test_data['reference']['msh'], checkIfExists: true) )
}