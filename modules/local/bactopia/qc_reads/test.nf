#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_READS } from './main.nf' 

workflow test_qc_reads_pe {

    inputs = tuple(
        [id:"output", runtype:"paired-end"],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC_READS ( inputs )
}

workflow test_qc_reads_se {

    inputs = tuple(
        [id:"output", runtype:"single-end"],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC_READS ( inputs )
}

workflow test_qc_reads_nanopore {

    inputs = tuple(
        [id:"output", runtype:"ont"],
        [file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC_READS ( inputs )
}
