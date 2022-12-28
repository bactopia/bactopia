#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STECFINDER } from './main.nf'
include { STECFINDER as STECFINDER_PE } from './main.nf' addParams(stecfinder_use_reads: true)

workflow test_stecfinder {

    inputs = tuple(
        [ id:"GCF_001695515", single_end:true, is_compressed:true ],
        file(params.test_data['species']['escherichia_coli']['genome']['fna_gz'], checkIfExists: true)
    )

    STECFINDER(inputs)
}

workflow test_stecfinder_pe {

    inputs = tuple(
        [ id:"ERR6005894", single_end:false, is_compressed:true, runtype:"illumina" ],
        [ file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]
    )

    STECFINDER_PE(inputs)
}
