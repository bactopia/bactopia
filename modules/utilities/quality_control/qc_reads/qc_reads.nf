
nextflow.enable.dsl = 2

process QC_READS {
    /* Cleanup the reads using Illumina-Cleanup */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*error.txt"

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra), path(genome_size)

    output:
    file "*-error.txt" optional true
    file "quality-control/*"
    tuple val(sample), val(single_end),
        path("quality-control/${sample}*.fastq.gz"),emit: READS,optional: true//,emit: COUNT_31MERS, ARIBA_ANALYSIS,MINMER_SKETCH, CALL_VARIANTS,MAPPING_QUERY optional true
    tuple val(sample), val(sample_type), val(single_end),
        path("quality-control/${sample}*.fastq.gz"), path(extra),
        path(genome_size),emit: ASSEMBLY, optional: true

    tuple val(sample), val(single_end),
        path("quality-control/${sample}*.{fastq,error-fq}.gz"),
        path(genome_size),emit: QC_FINAL_SUMMARY, optional: true
    file "${task.process}/*" optional true

    shell:
    qc_ram = task.memory.toString().split(' ')[0]
    is_assembly = sample_type.startsWith('assembly') ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapters = params.adapters ? path(params.adapters) : 'adapters'
    phix = params.phix ? path(params.phix) : 'phix'

    template "qc_reads.sh"

    stub:
    """
    mkdir quality-control
    mkdir ${task.process}
    touch ${sample}-error.txt
    touch quality-control/${sample}.fastq.gz
    touch quality-control/${sample}.error-fq.gz
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
        path(params.extra),
        path(params.genome_size)
    ])
    qc_reads(TEST_PARAMS_CH)
}
