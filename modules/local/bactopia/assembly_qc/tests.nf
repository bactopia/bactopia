#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ASSEMBLY_QC } from './main.nf' 

workflow test_assembly_qc {

    inputs = tuple(
        "test_assembly_qc",
        file(params.test_data['reference']['genome_size'], checkIfExists: true)
        file(params.test_data['reference']['fna'], checkIfExists: true),
        file(params.test_data['reference']['total_contigs', checkIfExists: true)
    )

    ASSEMBLY_QC ( inputs, ['checkm', 'quast'] )
}
