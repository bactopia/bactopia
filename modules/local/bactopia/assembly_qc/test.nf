#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ASSEMBLY_QC } from './main.nf' 

workflow test_assembly_qc {

    inputs = tuple(
        [id:"output"],
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['fna'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['total_contigs'], checkIfExists: true)
    )

    ASSEMBLY_QC ( inputs )
}
