#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SHIGATYPER } from './main.nf'

workflow test_shigatyper_pe {

    inputs = tuple(
        [ id:"ERR6005894", single_end:false, runtype:"illumina" ],
        [ file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]
    )

    SHIGATYPER ( inputs )
}

workflow test_shigatyper_se {

    inputs = tuple(
        [ id:"ERR6005894SE", single_end:true, runtype:"illumina" ],
        [ file(params.test_data['species']['shigella_dysenteriae']['illumina']['se'], checkIfExists: true) ]
    )

    SHIGATYPER ( inputs )
}

workflow test_shigatyper_ont {

    inputs = tuple(
        [ id:"SRR13039589", single_end:true, runtype:"ont" ],
        [ file(params.test_data['species']['shigella_dysenteriae']['nanopore']['se'], checkIfExists: true) ]
    )

    SHIGATYPER ( inputs )
}
