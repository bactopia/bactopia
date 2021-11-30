#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ISMAPPER } from './main.nf' 

workflow test_ismapper {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )
    reference = [file(params.test_data['species']['portiera']['genome']['gbk_gz'], checkIfExists: true)]
    query = [file(params.test_data['species']['haemophilus_influenzae']['genome']['is1016'], checkIfExists: true)]

    ISMAPPER ( inputs, reference, query )
}
