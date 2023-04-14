#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SKETCHER } from './main.nf' 

workflow test_sketcher_pe {
    inputs = tuple(
        [id:"output", single_end:false],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]
    )
    mash_db = file(params.test_data['datasets']['mash'], checkIfExists: true)
    sourmash_db = file(params.test_data['datasets']['sourmash'], checkIfExists: true)

    SKETCHER ( inputs, mash_db, sourmash_db )
}

workflow test_sketcher_se {
    inputs = tuple(
        [id:"output", single_end:true],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)],
        file(params.test_data['datasets']['mash'], checkIfExists: true),
        file(params.test_data['datasets']['sourmash'], checkIfExists: true)
    )
    mash_db = file(params.test_data['datasets']['mash'], checkIfExists: true)
    sourmash_db = file(params.test_data['datasets']['sourmash'], checkIfExists: true)

    SKETCHER ( inputs, mash_db, sourmash_db )
}
