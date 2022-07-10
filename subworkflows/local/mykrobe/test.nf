#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MYKROBE } from './main.nf' addParams(mykrobe_species: "staph")

workflow test_mykrobe {

    inputs = tuple(
        [ id:"test", runtype: "paired-end" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    MYKROBE(inputs)
}
