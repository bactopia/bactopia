#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ANNOTATE_GENOME } from './main.nf' 

workflow test_annotate_genome {

    inputs = tuple(
        [ id:"${params.test_data['reference']['name']}" ],
        file(params.test_data['reference']['genome_size']),
        file(params.test_data['reference']['fna_gz']),
        file(params.test_data['reference']['total_contigs'])
    )

    ANNOTATE_GENOME ( inputs, file(params.test_data['empty']['proteins']), file(params.test_data['empty']['prodigal_tf']) )
}

workflow test_annotate_genome_uncompressed {

    inputs = tuple(
        [ id:"${params.test_data['reference']['name']}" ],
        file(params.test_data['reference']['genome_size']),
        file(params.test_data['reference']['fna']),
        file(params.test_data['reference']['total_contigs'])
    )

    ANNOTATE_GENOME ( inputs, file(params.test_data['empty']['proteins']), file(params.test_data['empty']['prodigal_tf']) )
}
