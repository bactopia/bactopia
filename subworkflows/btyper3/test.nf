#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BTYPER3 } from './main.nf' 

workflow test_btyper3 {

    inputs = tuple(
        [ id:"GCF_000008445" ],
        file(params.test_data['species']['bacillus_anthracis']['genome']['fna_gz'], checkIfExists: true)
    )

    BTYPER3 ( inputs )
}
