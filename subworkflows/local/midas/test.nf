#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MIDAS } from './main.nf' addParams(midas_db: params.test_data['datasets']['midas'])
include { MIDAS as MIDAS_TARBALL } from './main.nf' addParams(midas_db: params.test_data['datasets']['midas_tarball'])

workflow test_midas {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    MIDAS( inputs )
}

workflow test_midas_tarball {

    inputs = tuple(
        [ id:"SRR2838702", single_end: true],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)]
    )

    MIDAS_TARBALL( inputs )
}
