#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATHER_SAMPLES } from './main.nf' 

workflow test_gather_samples_pe {

    inputs = tuple(
        'test_gather_samples_pe',
        'paired-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_pe {

    inputs = tuple(
        'test_gather_samples_merge_pe',
        'merge-pe',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_pe_hybrid {

    inputs = tuple(
        'test_gather_samples_merge_pe_hybrid',
        'hybrid-merge-pe',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)],
        file(params.test_data['nanopore']['se'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_se {

    inputs = tuple(
        'test_gather_samples_se',
        'single-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_merge_se {

    inputs = tuple(
        'test_gather_samples_merge_se',
        'merge-se',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_sra_accession {

    inputs = tuple(
        params.test_data['accessions']['srx'],
        'sra-accession',
        params.genome_size,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_assembly_accession {

    inputs = tuple(
        params.test_data['accessions']['gcf'],
        'assembly-accession',
        params.genome_size,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_assembly {

    inputs = tuple(
        params.test_data['reference']['name'],
        'assembly',
        params.genome_size,
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['reference']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}


workflow test_gather_samples_error_low_proportion {

    inputs = tuple(
        'test_gather_samples_error_low_proportion',
        'paired-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2_low_proportion'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_error_low_read_count_pe {

    inputs = tuple(
        'test_gather_samples_error_low_read_count_pe',
        'paired-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1_low_read_count'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2_low_read_count'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_error_low_read_count_se {

    inputs = tuple(
        'test_gather_samples_error_low_read_count_se',
        'single-end',
        params.genome_size,
        [file(params.test_data['illumina']['se_low_read_count'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_error_low_depth_pe {

    inputs = tuple(
        'test_gather_samples_error_low_depth_pe',
        'paired-end',
        params.genome_size,
        [file(params.test_data['illumina']['r1_low_depth'], checkIfExists: true)],
        [file(params.test_data['illumina']['r2_low_depth'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}

workflow test_gather_samples_error_low_depth_se {

    inputs = tuple(
        'test_gather_samples_error_low_depth_se',
        'single-end',
        params.genome_size,
        [file(params.test_data['illumina']['se_low_depth'], checkIfExists: true)],
        [file(params.test_data['empty']['fastq'], checkIfExists: true)],
        file(params.test_data['empty']['fna'], checkIfExists: true)
    )

    GATHER_SAMPLES ( inputs )
}
