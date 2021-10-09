#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MINMER_SKETCH } from './main.nf' 

workflow test_minmer_sketch_pe {
    inputs = tuple(
        [id:"output", single_end:false],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}

workflow test_minmer_sketch_se {
    inputs = tuple(
        [id:"output", single_end:true],
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}
