nextflow.enable.dsl = 2

process MAPPING_QUERY {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "mapping/*"

    input:
    tuple val(sample), val(single_end), path(fq)
    each path(query)

    output:
    file "mapping/*"
    file "${task.process}/*" optional true

    when:
    MAPPING_FASTAS.isEmpty() == false

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    template "mapping_query.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir mapping
    touch ${task.process}/${sample}
    touch mapping/${sample}
    """
}

//###############
//Module testing
//###############

workflow test{
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.query)
        )
    mapping_query(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect())
}
