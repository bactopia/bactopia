nextflow.enable.dsl = 2

process BLAST_GENES {
    /*
    Query gene FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "genes/*.{json,json.gz}"

    input:
    tuple val(sample), path(blastdb)
    path(query)

    output:
    path("genes/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_GENE_FASTAS.isEmpty() == false

    shell:
    template "blast_genes.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir genes
    touch ${task.process}/${sample}
    touch genes/${sample}.json
    touch genes/${sample}.json.gz
    """
}

//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        path(params.blastdb),
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.query)
        )

    blast_genes(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
