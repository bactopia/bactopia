#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ASSEMBLER } from './main.nf' addParams(reassemble: false, min_genome_size: 0)

workflow test_assembler_se {
    // [sample, runtype, single_end], fastqs, extra, genome_size
    inputs = tuple(
        [ id:"portiera", runtype:'single-end', single_end:true, genome_size: 350000 ],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ASSEMBLER ( inputs )
}

workflow test_assembler_pe {
    inputs = tuple(
        [ id:"portiera", runtype:'paired-end', single_end:false, genome_size:350000 ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ASSEMBLER ( inputs )
}

workflow test_assembler_hybrid {
    inputs = tuple(
        [ id:"portiera", runtype:'hybrid', single_end:false, genome_size:350000 ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true),
         file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ASSEMBLER ( inputs )
}

workflow test_assembler_short_polish {
    inputs = tuple(
        [ id:"portiera", runtype:'short_polish', single_end:true, genome_size:350000 ],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true),
         file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ASSEMBLER ( inputs )
}

workflow test_assembler_nanopore {
    inputs = tuple(
        [ id:"portiera", runtype:'ont', single_end:true, genome_size:350000 ],
        [file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna'])
    )

    ASSEMBLER ( inputs )
}
