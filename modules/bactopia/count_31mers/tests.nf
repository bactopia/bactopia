#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { COUNT_31MERS } from './main.nf' 

workflow test_count_31mers_pe {

    inputs = tuple(
        "test_count_31mers_pe",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    COUNT_31MERS ( inputs )
}

workflow test_count_31mers_se {

    inputs = tuple(
        "test_count_31mers_se",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    COUNT_31MERS ( inputs )
}
