#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MINMER_SKETCH } from './main.nf' 

workflow test_minmer_sketch_pe {
    inputs = tuple(
        "test_minmer_sketch_pe",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}

workflow test_minmer_sketch_se {
    inputs = tuple(
        "test_minmer_sketch_se",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}
