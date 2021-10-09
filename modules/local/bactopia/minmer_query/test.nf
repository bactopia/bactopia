#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MINMER_QUERY } from './main.nf' 

workflow test_mash_query_pe {
    inputs = tuple(
        [ id:"output", single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['illumina']['sig'], checkIfExists: true)
    )

    MINMER_QUERY ( inputs, Channel.fromPath("${params.test_data['datasets']['minmer']}/*.msh") )
}

workflow test_mash_query_se {
    inputs = tuple(
        [ id:"output", single_end:true ],
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['illumina']['sig'], checkIfExists: true)
    )

    MINMER_QUERY ( inputs, Channel.fromPath("${params.test_data['datasets']['minmer']}/*.msh") )
}

workflow test_sourmash_query_pe {
    inputs = tuple(
        [ id:"output", single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['illumina']['sig'], checkIfExists: true)
    )

    MINMER_QUERY ( inputs, Channel.fromPath("${params.test_data['datasets']['minmer']}/*.json.gz") )
}

workflow test_sourmash_query_se {
    inputs = tuple(
        [ id:"output", single_end:true ],
        [file(params.test_data['illumina']['se'], checkIfExists: true)],
        file(params.test_data['illumina']['sig'], checkIfExists: true)
    )

    MINMER_QUERY ( inputs, Channel.fromPath("${params.test_data['datasets']['minmer']}/*.json.gz") )
}

