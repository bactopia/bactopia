#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SNIPPY } from './main.nf' addParams(reference: params.test_data['species']['portiera']['genome']['gbk_gz'])
include { SNIPPY as SNIPPY_MASK } from './main.nf' addParams(reference: params.test_data['species']['portiera']['genome']['gbk_gz'],
                                                             mask: params.test_data['species']['portiera']['genome']['bed'])

workflow test_snippy {
    inputs = Channel.of(
        [[ id:"test0", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test1", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test2", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test3", runtype: "paired-end" ],
         [file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test4", runtype: "paired-end" ],
         [file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]]
    )

    SNIPPY(inputs)
}

workflow test_snippy_mask {
    inputs = Channel.of(
        [[ id:"test0", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test1", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test2", runtype: "paired-end" ],
         [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test3", runtype: "paired-end" ],
         [file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]],
        [[ id:"test4", runtype: "paired-end" ],
         [file(params.test_data['species']['shigella_dysenteriae']['illumina']['r1'], checkIfExists: true),
          file(params.test_data['species']['shigella_dysenteriae']['illumina']['r2'], checkIfExists: true)]]
    )

    SNIPPY_MASK(inputs)
}
