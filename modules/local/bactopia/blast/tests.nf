#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLAST } from './main.nf' 

workflow test_blast_genes {

    inputs = tuple(
        [id:params.test_data['reference']['name']],
        file("${params.test_data['reference']['blastdb']}/*")
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['genes'] )
}

workflow test_blast_primers {

    inputs = tuple(
        [id:params.test_data['reference']['name']],
        file("${params.test_data['reference']['blastdb']}/*")
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['primers'] )
}

workflow test_blast_proteins {

    inputs = tuple(
        [id:params.test_data['reference']['name']],
        file("${params.test_data['reference']['blastdb']}/*")
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['proteins'] )
}
