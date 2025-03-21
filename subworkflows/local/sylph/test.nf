#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SYLPH } from './main.nf' addParams(sylph_db: params.test_data['datasets']['sylph'])

workflow test_sylph {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    SYLPH ( inputs )
} 
