#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MAPPING_QUERY } from './main.nf' 

workflow test_mapping_query_pe {

    inputs = tuple(
        "test_mapping_query_pe",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    MAPPING_QUERY ( inputs, params.test_data['datasets']['mapping'] )
}

workflow test_mapping_query_se {
    inputs = tuple(
        "test_mapping_query_se",
        false,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    MAPPING_QUERY ( inputs, params.test_data['datasets']['mapping']  )
}