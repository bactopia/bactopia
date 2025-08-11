#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { AMRFINDERPLUS } from './main.nf' 

workflow test_amrfinderplus {

    inputs = tuple(
        [ id:"SRR2838702" ],
        file(params.test_data['species']['portiera']['illumina']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['illumina']['faa_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['illumina']['gff_gz'], checkIfExists: true),
    )

    db = file(params.test_data['datasets']['amrdb']['amrfinder'], checkIfExists: true)

    AMRFINDERPLUS ( inputs, db )
}
