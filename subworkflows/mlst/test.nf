#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MLST } from './main.nf' 

workflow test_mlst {

    inputs = tuple(
        [ id:"GCF_000017085" ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna'], checkIfExists: true)
    )

    db = file(params.test_data['datasets']['mlst'], checkIfExists: true)

    MLST ( inputs, db )
}
