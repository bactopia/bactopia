#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EMMTYPER } from './main.nf'
include { EMMTYPER as EMMTYPER_WITH_DB } from './main.nf' addParams(emmtyper_blastdb: params.test_data['datasets']['emmtyper'])

workflow test_emmtyper {

    inputs = tuple(
        [ id:"GCF_006364235" ],
        file(params.test_data['species']['streptococcus_pyogenes']['genome']['fna_gz'], checkIfExists: true)
    )

    EMMTYPER ( inputs )
}

workflow test_emmtyper_db {

    inputs = tuple(
        [ id:"GCF_006364235" ],
        file(params.test_data['species']['streptococcus_pyogenes']['genome']['fna_gz'], checkIfExists: true)
    )

    EMMTYPER_WITH_DB ( inputs )
}
