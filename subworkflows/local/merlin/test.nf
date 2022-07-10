#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MERLIN } from './main.nf' addParams(mash_sketch: params.test_data['datasets']['refseq_msh'])

workflow test_merlin {

    inputs = tuple(
        [ id:"GCF_000017085", single_end: false ],
        file(params.test_data['species']['staphylococcus_aureus']['genome']['fna_gz'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    MERLIN(inputs)
}

workflow test_merlinfull {

    inputs = tuple(
        [ id:"GCF_000292685", single_end: false ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    MERLIN(inputs)
}
