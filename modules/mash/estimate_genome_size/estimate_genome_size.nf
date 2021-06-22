nextflow.enable.dsl = 2

process ESTIMATE_GENOME_SIZE {
    /* Estimate the input genome size if not given. */
    tag "${sample}"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra)

    output:
    path "${sample}-genome-size-error.txt" optional true
    path("${sample}-genome-size.txt") optional true
    tuple val(sample), val(sample_type), val(single_end), 
        path("fastqs/${sample}*.fastq.gz"), path(extra), path("${sample}-genome-size.txt"),emit: QUALITY_CONTROL, optional: true
    path "${task.process}/*" optional true

    shell:
    genome_size = SPECIES_GENOME_SIZE
    
    template "estimate_genome_size.sh"

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch ${sample}-genome-size-error.txt
    touch ${sample}-genome-size.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        path(params.fq),
        path(params.extra)             
        ])

    estimate_genome_size(TEST_PARAMS_CH)
}
