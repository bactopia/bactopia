#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQUENCE_TYPE } from './main.nf' 

workflow test_sequence_type_pe {

    inputs = tuple(
        [id:"test_sequence_type_pe", single_end:false],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['reference']['fna_gz'], checkIfExists: true)
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_pe_uncompressed {

    inputs = tuple(
        [id:"test_sequence_type_pe", single_end:false],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_se {

    inputs = tuple(
        [id:"test_sequence_type_se", single_end:true],
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['reference']['fna_gz'], checkIfExists: true)
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_se_uncompressed {

    inputs = tuple(
        [id:"test_sequence_type_se", single_end:true],
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}
