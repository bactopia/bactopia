#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCOARY } from './main.nf' addParams(traits: params.test_data['species']['portiera']['genome']['gene_traits'])

workflow test_scoary {

    inputs = tuple(
        [ id:"scoary" ],
        file(params.test_data['species']['portiera']['genome']['gene_presence_absence'], checkIfExists: true)
    )

    SCOARY ( inputs )
}
