nextflow.enable.dsl = 2

process GATHER_FASTQS {
    /* Gather up input FASTQs for analysis. */
    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "bactopia.versions"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    tag "${sample}"

    input:
    tuple val(sample), val(sample_type), val(single_end), path(r1: '*???-r1'), path(r2: '*???-r2'), path(extra)

    output:
    path("*-error.txt") optional true
    tuple val(sample), val(final_sample_type), val(single_end),
        path("fastqs/${sample}*.fastq.gz"), path("extra/*.gz"), emit: FASTQ_PE_STATUS, optional: true
    path("${task.process}/*") optional true
    path("bactopia.versions") optional true
    path("multiple-read-sets-merged.txt") optional true

    shell:
    bactopia_version = VERSION
    nextflow_version = nextflow.version
    is_assembly = sample_type.startsWith('assembly') ? true : false
    is_compressed = false
    no_cache = params.no_cache ? '-N' : ''
    use_ena = params.use_ena
    if (task.attempt >= 4) {
        if (use_ena) {
            // Try SRA
            use_ena = false 
        } else {
            // Try ENA
            use_ena = true
        }
    }
    if (extra) {
        is_compressed = extra.getName().endsWith('gz') ? true : false
    }
    section = null
    if (sample_type == 'assembly_accession') {
        section = sample.startsWith('GCF') ? 'refseq' : 'genbank'
    }
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    final_sample_type = sample_type
    if (sample_type == 'hybrid-merge-pe') {
        final_sample_type = 'hybrid'
    } else if (sample_type == 'merge-pe') {
        final_sample_type = 'paired-end'
    } else if (sample_type == 'merge-se') {
        final_sample_type = 'single-end'
    }

    template "gather_fastqs.sh"

    stub:
    final_sample_type = 'single-end'
    """
    mkdir fastqs
    mkdir extra
    mkdir ${task.process}
    touch ${sample}-error.txt
    touch fastqs/${sample}.fastq.gz
    touch extra/${sample}.gz
    touch ${task.process}/${sample}
    touch bactopia.versions
    touch multiple-read-sets-merged.txt
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    test_params_input = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        params.r1,
        params.r2,
        params.extra           
        ])

    gather_fastqs(test_params_input)
}
