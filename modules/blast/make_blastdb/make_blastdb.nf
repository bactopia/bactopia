nextflow.enable.dsl = 2

process MAKE_BLASTDB {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    tuple val(sample), val(single_end), path(fasta)

    output:
    path("blastdb/*")
    tuple val(sample), path("blastdb/*"), emit: BLAST_DB, optional:true
    file "${task.process}/*" optional true

    shell:
    template "make_blastdb.sh"

    stub:
    """
    mkdir blastdb
    mkdir ${task.process}
    touch blastdb/${sample}
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
        path(params.fasta)
    ])

    make_blastdb(TEST_PARAMS_CH)
}
