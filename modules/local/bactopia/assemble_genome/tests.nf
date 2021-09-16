#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ASSEMBLE_GENOME } from './main.nf' 

workflow test_assemble_genome_se {
    // [sample, runtype, single_end], fastqs, extra, genome_size
    inputs = tuple(
        [ id:"test_assemble_genome_se", runtype:'single-end', single_end:true ],
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['empty']['fna']),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_pe {
    // [sample, runtype, single_end], fastqs, extra, genome_size
    inputs = tuple(
        [ id:"test_assemble_genome_pe", runtype:'paired-end', single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna']),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}

workflow test_assemble_genome_hybrid {
    // [sample, runtype, single_end], fastqs, extra, genome_size
    inputs = tuple(
        [ id:"test_assemble_genome_hybrid", runtype:'hybrid', single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['nanopore']['se'], checkIfExists: true),
        file(params.test_data['reference']['genome_size'])
    )

    ASSEMBLE_GENOME ( inputs )
}
