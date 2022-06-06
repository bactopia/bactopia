#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ARIBA } from './main.nf' 

workflow test_ariba {

    inputs = tuple(
        [ id:"SRR2838702", single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    ARIBA ( inputs )
}
