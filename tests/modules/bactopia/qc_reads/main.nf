#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANNOTATE_GENOME } from '../../../modules/bactopia/annotate_genome/main.nf' 

workflow test_annotate_genome {

    inputs = tuple(
        sample,
        false,
        [file(params.test_data[''][''][''], checkIfExists: true), file(params.test_data[''][''][''], checkIfExists: true)],
        file(fasta),
        file(total_contigs)
    )

    ANNOTATE_GENOME ( inputs, file(prokka_proteins), file(prodigal_tf) )
}
