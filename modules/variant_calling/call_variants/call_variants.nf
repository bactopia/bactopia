nextflow.enable.dsl = 2

process CALL_VARIANTS {
    /*
    Identify variants (SNPs/InDels) against a set of reference genomes
    using Snippy.
    */
    tag "${sample} - ${reference_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/user", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    tuple val(sample), val(single_end), path(fq)
    each path(reference)

    output:
    path "${reference_name}/*"
    path "${task.process}/*" optional true

    when:
    REFERENCES.isEmpty() == false

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template "call_variants.sh"

    stub:
    reference_name = reference.getSimpleName()
    """
    mkdir ${reference_name}
    mkdir ${task.process}
    touch ${reference_name}/*
    touch ${task.process}/*
    """
}

//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq),
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.reference)
        )
    call_variants(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect)
}
