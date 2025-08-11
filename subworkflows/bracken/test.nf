#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BRACKEN } from './main.nf' addParams(kraken2_db: params.test_data['datasets']['kraken2'])
include { BRACKEN as BRACKEN_TARBALL } from './main.nf' addParams(kraken2_db: params.test_data['datasets']['kraken2_tarball'])

workflow test_bracken {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    BRACKEN( inputs )
}

workflow test_bracken_tarball {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    BRACKEN_TARBALL( inputs )
}
