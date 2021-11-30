#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EMMTYPER } from './main.nf' 

workflow test_emmtyper {

    inputs = tuple(
        [ id:"GCF_006364235" ],
        file(params.test_data['species']['streptococcus_pyogenes']['genome']['fna_gz'], checkIfExists: true)
    )

    EMMTYPER ( inputs )
}
