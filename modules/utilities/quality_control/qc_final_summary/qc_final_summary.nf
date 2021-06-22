nextflow.enable.dsl = 2

process QC_FINAL_SUMMARY {
    /* Run FASTQC on the cleaned up FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    tuple val(sample), val(single_end), path(fq), path(genome_size)

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:

    template "qc_final_summary.sh"

    stub:
    """
    mkdir quality-control
    mkdir ${task.process}
    touch quality-control/${sample}
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing
//###############

workflow test{

    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq),
        path(params.genome_size)
        ])

    qc_final_summary(TEST_PARAMS_CH)
}
