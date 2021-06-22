nextflow.enable.dsl = 2

process DOWNLOAD_REFERENCES {
    /*
    Download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${params.max_references} reference(s)"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: 'mash-dist.txt'

    input:
    tuple val(sample), val(single_end), path(fq), path(sample_sketch)
    path(refseq_sketch)

    output:
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("genbank/*.gbk"), emit:CALL_VARIANTS_AUTO, optional: true
    path("mash-dist.txt")
    file "${task.process}/*" optional true

    when:
    REFSEQ_SKETCH_FOUND == true

    shell:
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references
    template "download_references.sh"

    stub:
    """
    mkdir fastqs
    mkdir genbank
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch genbank/*.gbk
    touch ${task.process}/${sample}
    touch mash-dist.txt
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
        path(params.sample_sketch)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.refseq_sketch)
        )
    download_references(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}

