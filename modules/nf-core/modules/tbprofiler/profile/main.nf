// Import generic module functions
include { initOptions; saveFiles } from '../../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options, 'tbprofiler')
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::tb-profiler=3.0.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:3.0.8--pypyh5e36f6f_0' :
        'quay.io/biocontainers/tb-profiler:3.0.8--pypyh5e36f6f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bam/*.bam")      , emit: bam
    tuple val(meta), path("results/*.csv")  , emit: csv, optional: true
    tuple val(meta), path("results/*.json") , emit: json
    tuple val(meta), path("results/*.txt")  , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")   , emit: vcf
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    """
    tb-profiler \\
        profile \\
        $options.args \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        $input_reads

    cat <<-END_VERSIONS > versions.yml
    tbprofiler:
        tb-profiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}
