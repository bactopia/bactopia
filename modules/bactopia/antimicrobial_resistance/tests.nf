#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANTIMICROBIAL_RESISTANCE } from './main.nf'

workflow test_antimicrobial_resistance {

    // sample, genes, proteins (e.g. GCF_000292685, GCF_000292685.ffn, GCF_000292685.faa)
    inputs = tuple(
        params.test_data['reference'],
        file(params.test_data['reference_ffn'], checkIfExists: true),
        file(params.test_data['reference_faa'], checkIfExists: true)
    )

    amrdbs = [
        file(params.test_data['databases']['amrdb']['amrfinder'], checkIfExists: true)
    ]

    ANTIMICROBIAL_RESISTANCE ( inputs, amrdbs )
}
