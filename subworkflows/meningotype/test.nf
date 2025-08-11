#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MENINGOTYPE } from './main.nf' 

workflow test_meningotype {

    inputs = tuple(
        [ id:"GCF_003355215" ],
        file(params.test_data['species']['neisseria_meningitidis']['genome']['fna_gz'], checkIfExists: true)
    )

    MENINGOTYPE ( inputs )
}
