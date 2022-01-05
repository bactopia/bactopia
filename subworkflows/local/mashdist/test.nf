#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASHDIST } from './main.nf' addParams(mash_sketch: params.test_data['datasets']['refseq_genomes'])

workflow test_mashdist {

    inputs = tuple(
        [ id:"GCF_000292685" ],
        file(params.test_data['species']['portiera']['genome']['fna_gz'], checkIfExists: true)
    )

    MASHDIST(inputs)
}
