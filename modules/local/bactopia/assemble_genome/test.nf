#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ASSEMBLE_GENOME } from './main.nf' 

workflow test_assemble_genome_se {
    // [sample, runtype, single_end], fastqs, extra, genome_size
    inputs = tuple(
        [ id:"output", runtype:'single-end', single_end:true ],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
        file(params.test_data['empty']['fna']),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_pe {
    inputs = tuple(
        [ id:"output", runtype:'paired-end', single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
        file(params.test_data['empty']['fna']),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_hybrid {
    inputs = tuple(
        [ id:"output", runtype:'hybrid', single_end:false ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
        file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_nanopore {
    inputs = tuple(
        [ id:"output", runtype:'ont', single_end:true ],
        [file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna']),
        file(params.test_data['species']['portiera']['genome']['genome_size'], checkIfExists: true)
    )

    ASSEMBLE_GENOME ( inputs )
}
