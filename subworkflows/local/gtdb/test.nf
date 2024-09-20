#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GTDB } from './main.nf' addParams(gtdb: params.test_data['datasets']['gtdb'])
include { GTDB as GTDB_TARBALL } from './main.nf' addParams(gtdb: params.test_data['datasets']['gtdb_tarball'])

workflow test_gtdb {

    inputs = tuple(
        [ id:"test" ],
        [ file(params.test_data['species']['staphylococcus_aureus']['genome']['fna_gz'], checkIfExists: true),
          file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true) ]
    )

    GTDB ( inputs )
}

workflow test_gtdb_tarball {

    inputs = tuple(
        [ id:"test" ],
        [ file(params.test_data['species']['staphylococcus_aureus']['genome']['fna_gz'], checkIfExists: true),
          file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true) ]
    )

    GTDB_TARBALL ( inputs )
}
 