#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ARIBA_ANALYSIS } from './main.nf' 

workflow test_ariba_analysis {

    // sample, single_end, fastqs
    inputs = tuple(
        [ id:'output', single_end:false ],
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    ariba_dbs = [
        file(params.test_data['datasets']['ariba']['card'], checkIfExists: true),
        file(params.test_data['datasets']['ariba']['vfdb_core'], checkIfExists: true)
    ]
    ARIBA_ANALYSIS ( inputs, ariba_dbs )
}
