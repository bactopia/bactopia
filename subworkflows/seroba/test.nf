#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEROBA } from './main.nf' 

workflow test_seroba {

    inputs = tuple(
        [ id:"test", runtype: "paired-end", single_end: false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    SEROBA ( inputs )
}
