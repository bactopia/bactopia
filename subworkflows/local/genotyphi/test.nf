#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOTYPHI } from './main.nf' 

workflow test_genotyphi {

    inputs = tuple(
        [ id:"test", runtype: "paired-end" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    GENOTYPHI(inputs)
}
