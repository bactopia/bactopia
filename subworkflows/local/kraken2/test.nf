#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KRAKEN2 } from './main.nf' addParams(kraken2_db: params.test_data['datasets']['kraken2'])

workflow test_kraken2 {

    inputs = tuple(
        [ id:"SRR2838702" ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    KRAKEN2( inputs )
}
