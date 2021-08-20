#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CALL_VARIANTS } from './main.nf' 

workflow test_call_variants_pe {

    inputs = tuple(
        params.test_data['illumina']['name'],
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    references = [
        file(params.test_data['reference']['gbk'], checkIfExists: true)
    ]

    CALL_VARIANTS ( inputs, references )
}

workflow test_call_variants_se {

    inputs = tuple(
        params.test_data['illumina']['name'],
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    references = [
        file(params.test_data['reference']['gbk'], checkIfExists: true)
    ]

    CALL_VARIANTS ( inputs, references )
}

workflow test_call_variants_auto_pe {

    inputs = tuple(
        params.test_data['illumina']['name'],
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    references = [
        file(params.test_data['datasets']['refseq_genomes'], checkIfExists: true)
    ]

    CALL_VARIANTS ( inputs, references )
}

workflow test_call_variants_auto_se {

    inputs = tuple(
        params.test_data['illumina']['name'],
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    references = [
        file(params.test_data['datasets']['refseq_genomes'], checkIfExists: true)
    ]

    CALL_VARIANTS ( inputs, references )
}
