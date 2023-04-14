#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATHER } from './main.nf' 

workflow test_gather_pe {

    inputs = tuple(
        [id:'output', runtype:'paired-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )
    GATHER ( inputs )
}

workflow test_gather_merge_pe {

    inputs = tuple(
        [id:'output', runtype:'merge-pe', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_merge_pe_hybrid {

    inputs = tuple(
        [id:'output', runtype:'hybrid-merge-pe', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['species']['portiera']['nanopore']['se'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_se {

    inputs = tuple(
        [id:'output', runtype:'single-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['se'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_merge_se {

    inputs = tuple(
        [id:'output', runtype:'merge-se', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true),
         file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_sra_accession {

    inputs = tuple(
        [id:params.test_data['accessions']['srx'], runtype:'sra_accession', genome_size:params.genome_size],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_assembly_accession {

    inputs = tuple(
        [id:params.test_data['accessions']['gcf'], runtype:'assembly_accession', genome_size:params.genome_size],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_assembly {

    inputs = tuple(
        [id:params.test_data['species']['portiera']['genome']['name'], runtype:'assembly', genome_size:params.genome_size],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['species']['portiera']['genome']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}


workflow test_gather_error_low_proportion {

    inputs = tuple(
        [id:'output', runtype:'paired-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2_low_proportion'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_error_low_read_count_pe {

    inputs = tuple(
        [id:'output', runtype:'paired-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1_low_read_count'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2_low_read_count'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_error_low_read_count_se {

    inputs = tuple(
        [id:'output', runtype:'single-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['se_low_read_count'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_error_low_depth_pe {

    inputs = tuple(
        [id:'output', runtype:'paired-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['r1_low_depth'], checkIfExists: true)],
        [file(params.test_data['species']['portiera']['illumina']['r2_low_depth'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}

workflow test_gather_error_low_depth_se {

    inputs = tuple(
        [id:'output', runtype:'single-end', genome_size:params.genome_size],
        [file(params.test_data['species']['portiera']['illumina']['se_low_depth'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER ( inputs )
}
