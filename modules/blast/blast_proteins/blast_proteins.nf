nextflow.enable.dsl = 2

process BLAST_PROTEINS {
    /*
    Query protein FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "proteins/*.{json,json.gz}"

    input:
    tuple val(sample), path(blastdb)
    path(query)

    output:
    path("proteins/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PROTEIN_FASTAS.isEmpty() == false

    shell:

    template "blast_proteins.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir proteins
    touch ${task.process}/${sample}
    touch proteins/${sample}.json
    touch proteins/${sample}.json.gz
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

    blast_proteins(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
