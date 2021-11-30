#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANNOTATE_GENOME } from './main.nf' 

workflow test_annotate_genome {

    inputs = tuple(
        [ id:"${params.test_data['species']['portiera']['genome']['name']}" ],
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['total_contigs'], checkIfExists: true)
    )

    ANNOTATE_GENOME ( inputs, file(params.test_data['empty']['proteins']), file(params.test_data['empty']['prodigal_tf']) )
}

workflow test_annotate_genome_uncompressed {

    inputs = tuple(
        [ id:"${params.test_data['species']['portiera']['genome']['name']}" ],
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['fna'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['total_contigs'], checkIfExists: true)
    )

    ANNOTATE_GENOME ( inputs, file(params.test_data['empty']['proteins']), file(params.test_data['empty']['prodigal_tf']) )
}
