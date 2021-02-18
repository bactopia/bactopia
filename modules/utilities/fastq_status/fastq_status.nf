nextflow.enable.dsl = 2

process FASTQ_STATUS {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra)
    output:
    file "*-error.txt" optional true
    tuple val(sample), val(sample_type), val(single_end), 
        path("fastqs/${sample}*.fastq.gz"), path(extra),emit: ESTIMATE_GENOME_SIZE, optional: true
    file "${task.process}/*" optional true

    shell:
    single_end = fq[1] == null ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    
    template "fastq_status.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir fastqs
    touch ${sample}-error.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        path(params.fq),
        path(params.extra)             
        ])

    fastq_status(TEST_PARAMS_CH)
}
