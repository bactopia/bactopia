#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANTIMICROBIAL_RESISTANCE } from './main.nf'

workflow test_antimicrobial_resistance {

    // sample, genes, proteins (e.g. GCF_000292685, GCF_000292685.ffn, GCF_000292685.faa)
    inputs = tuple(
        [ id:"${params.test_data['reference']['name']}" ],
        file(params.test_data['species']['portiera']['genome']['ffn_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['faa_gz'], checkIfExists: true)
    )

    amrdbs = [
        file(params.test_data['datasets']['amrdb']['amrfinder'], checkIfExists: true)
    ]

    ANTIMICROBIAL_RESISTANCE ( inputs, amrdbs )
}

workflow test_antimicrobial_resistance_uncompressed {

    // sample, genes, proteins (e.g. GCF_000292685, GCF_000292685.ffn, GCF_000292685.faa)
    inputs = tuple(
        [ id:"${params.test_data['reference']['name']}" ],
        file(params.test_data['species']['portiera']['genome']['ffn'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['faa'], checkIfExists: true)
    )

    amrdbs = [
        file(params.test_data['datasets']['amrdb']['amrfinder'], checkIfExists: true)
    ]

    ANTIMICROBIAL_RESISTANCE ( inputs, amrdbs )
}
