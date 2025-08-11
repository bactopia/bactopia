#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PHISPY } from './main.nf' 

workflow test_phispy {

    inputs = tuple(
        [ id:"GCF_900478275l" ],
        file(params.test_data['species']['haemophilus_influenzae']['genome']['gbk_gz'], checkIfExists: true)
    )

    PHISPY ( inputs )
}

