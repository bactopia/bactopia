#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLAST; MAKE_BLASTDB } from './main.nf' 

workflow test_blast_genes {

    inputs = tuple(
        params.test_data['reference']['name'],
        params.test_data['reference']['blastdb']
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['genes'] )
}

workflow test_blast_primers {

    inputs = tuple(
        params.test_data['reference']['name'],
        params.test_data['reference']['blastdb']
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['primers'] )
}

workflow test_blast_proteins {

    inputs = tuple(
        params.test_data['reference']['name'],
        params.test_data['reference']['blastdb']
    )

    BLAST ( inputs, params.test_data['datasets']['blast']['proteins'] )
}

workflow test_makeblastdb {

    inputs = tuple(
        params.test_data['reference']['name'],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    MAKE_BLASTDB ( inputs )
}
