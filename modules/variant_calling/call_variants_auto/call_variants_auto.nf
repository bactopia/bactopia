nextflow.enable.dsl = 2

process CALL_VARIANTS_AUTO {
    /*
    Identify variants (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample} - ${reference_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    tuple val(sample), val(single_end), path(fq), path(reference)

    output:
    path "${reference_name}/*"
    path "${task.process}/*" optional true

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName().split("${sample}-")[1].split(/\./)[0]
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template "call_variants_auto.sh"

    stub:
    reference_name = "ref_name"
    """
    echo True
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
        path(params.reference)
        ])
    call_variants_auto(TEST_PARAMS_CH)
}
