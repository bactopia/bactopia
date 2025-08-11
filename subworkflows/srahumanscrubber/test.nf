#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { SRAHUMANSCRUBBER } from './main.nf' 

workflow test_srahumanscrubber_pe {

    inputs = tuple(
        [ id:"SRR2838702", single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    SRAHUMANSCRUBBER ( inputs )
}

workflow test_srahumanscrubber_se {

    inputs = tuple(
        [ id:"SRR2838702", single_end:true ],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
    )

    SRAHUMANSCRUBBER ( inputs )
}
