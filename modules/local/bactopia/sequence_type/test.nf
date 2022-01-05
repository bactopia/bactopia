#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQUENCE_TYPE } from './main.nf' 

workflow test_sequence_type_pe {

    inputs = tuple(
        [id:"output", single_end:false],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_pe_uncompressed {

    inputs = tuple(
        [id:"output", single_end:false],
        file(params.test_data['species']['portiera']['genome']['fna'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_se {

    inputs = tuple(
        [id:"output", single_end:true],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}

workflow test_sequence_type_se_uncompressed {

    inputs = tuple(
        [id:"output", single_end:true],
        file(params.test_data['species']['portiera']['genome']['fna'], checkIfExists: true),
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)]
    )

    SEQUENCE_TYPE ( inputs, Channel.fromPath("${params.test_data['datasets']['mlst']}/*.tar.gz") )
}
