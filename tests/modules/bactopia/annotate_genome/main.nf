#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANNOTATE_GENOME } from '../../../../modules/bactopia/annotate_genome/main.nf' 

workflow test_annotate_genome {

    inputs = tuple(
        "test_annotate_genome",
        false,
        [file(params.test_data['illumina_r1'], checkIfExists: true), file(params.test_data['illumina_r2'], checkIfExists: true)],
        file(params.test_data['reference_fna']),
        file(params.test_data['total_contigs'])
    )

    ANNOTATE_GENOME ( inputs, file(params.test_data['prokka_proteins']), file(params.test_data['prodigal_tf']) )
}
