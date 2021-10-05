#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MAPPING_QUERY } from './main.nf' 

workflow test_mapping_query_pe {

    inputs = tuple(
        [ id:"test_mapping_query_pe", single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    MAPPING_QUERY ( inputs, params.test_data['datasets']['mapping'] )
}

workflow test_mapping_query_se {
    inputs = tuple(
        [ id:"test_mapping_query_se", single_end:true ],
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    MAPPING_QUERY ( inputs, params.test_data['datasets']['mapping']  )
}
