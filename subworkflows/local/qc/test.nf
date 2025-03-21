#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC } from './main.nf' 

workflow test_qc_pe {

    inputs = tuple(
        [id:"output", runtype:"paired-end", genome_size:350000],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC ( inputs )
}

workflow test_qc_se {

    inputs = tuple(
        [id:"output", runtype:"single-end", genome_size:350000],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC ( inputs )
}

workflow test_qc_nanopore {

    inputs = tuple(
        [id:"output", runtype:"ont", genome_size:350000],
        [file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true),
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC ( inputs )
}

workflow test_qc_hybrid {

    inputs = tuple(
        [id:"output", runtype:"hybrid", genome_size:350000],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['adapters'], checkIfExists: true),
        file(params.test_data['empty']['phix'], checkIfExists: true)
    )

    QC ( inputs )
}
