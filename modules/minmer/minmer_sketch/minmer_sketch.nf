nextflow.enable.dsl = 2

process MINMER_SKETCH {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    tuple val(sample), val(single_end), path(fq)

    output:
    path("${sample}*.{msh,sig}")
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("${sample}.sig"),emit: MINMER_QUERY
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("${sample}-k31.msh"),emit: DOWNLOAD_REFERENCES
    path "${task.process}/*" optional true

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template "minmer_sketch.sh"

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/${sample}
    touch ${sample}.sig
    touch ${sample}-k31.msh

    """
}

//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq)
        ])

    minmer_sketch(TEST_PARAMS_CH)
}
