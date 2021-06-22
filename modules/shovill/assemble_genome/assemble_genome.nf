nextflow.enable.dsl = 2

process ASSEMBLE_GENOME {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "assembly/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${sample}-assembly-error.txt"

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra), path(genome_size)

    output:
    path "assembly/*"
    path "${sample}-assembly-error.txt" optional true
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("assembly/${sample}.{fna,fna.gz}"),emit: SEQUENCE_TYPE, optional:true
    tuple val(sample), val(single_end), path("assembly/${sample}.{fna,fna.gz}"), emit: MAKE_BLASTDB, optional: true
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("assembly/${sample}.{fna,fna.gz}"), path("total_contigs_*"),emit: ANNOTATION, optional:true
    tuple val(sample), path("assembly/${sample}.{fna,fna.gz}"), path(genome_size),emit: ASSEMBLY_QC, optional: true
    path "${task.process}/*" optional true

    shell:
    shovill_ram = task.memory.toString().split(' ')[0]
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    no_miniasm = params.no_miniasm ? "--no_miniasm" : ""
    no_rotate = params.no_rotate ? "--no_rotate" : ""
    no_pilon = params.no_pilon ? "--no_pilon" : ""
    keep = params.keep_all_files ? "--keep 3" : "--keep 1"
    use_original_assembly = null
    if (sample_type.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    template "assemble_genome.sh"

    stub:
    """
    mkdir assembly
    mkdir fastqs
    mkdir ${task.process}
    touch total_contigs_${sample}
    touch ${sample}-assembly-error.txt
    touch fastqs/${sample}.fastq.gz
    touch assembly/${sample}
    touch assembly/${sample}.fna
    touch assembly/${sample}.fna.gz
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

    assemble_genome(TEST_PARAMS_CH)
}
