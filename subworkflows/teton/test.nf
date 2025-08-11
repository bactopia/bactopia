#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TETON } from './main.nf' addParams(kraken2_db: params.test_data['datasets']['kraken2'])

workflow test_teton_pe {

    inputs = tuple(
        [ id:"SRR2838702", single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    TETON ( inputs )
}

workflow test_teton_se {

    inputs = tuple(
        [ id:"SRR2838702", single_end:true ],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
    )

    TETON ( inputs )
}
