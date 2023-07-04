#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PNEUMOCAT } from './main.nf' 

workflow test_pneumocat {

    inputs = tuple(
        [ id:"ERR1438863", runtype: "paired-end", single_end: false ],
        [file(params.test_data['species']['streptococcus_pneumoniae']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['streptococcus_pneumoniae']['illumina']['r2'], checkIfExists: true)]
    )

    PNEUMOCAT ( inputs )
}
