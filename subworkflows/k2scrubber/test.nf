#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { K2SCRUBBER } from './main.nf' 


include { SCRUBBER } from './main.nf' addParams(scrubber_db: params.test_data['datasets']['scrubber'])

workflow test_scrubber_pe {

    inputs = tuple(
        [ id:"SRR2838702", single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    K2SCRUBBER ( inputs )
}

workflow test_scrubber_se {

    inputs = tuple(
        [ id:"SRR2838702", single_end:true ],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
    )

    K2SCRUBBER ( inputs )
}
