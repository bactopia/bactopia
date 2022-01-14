// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'tbprofiler')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::tb-profiler=3.0.8"
conda_env   = file("${params.condadir}/tbprofiler").exists() ? "${params.condadir}/tbprofiler" : conda_tools

process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:3.0.8--pypyh5e36f6f_0' :
        'quay.io/biocontainers/tb-profiler:3.0.8--pypyh5e36f6f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bam/*.bam")     , emit: bam
    tuple val(meta), path("results/*.csv") , emit: csv, optional: true
    tuple val(meta), path("results/*.json"), emit: json
    tuple val(meta), path("results/*.txt") , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")  , emit: vcf
    path "*.{log,err}"                     , emit: logs, optional: true
    path ".command.*"                      , emit: nf_logs
    path "versions.yml"                    , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    def platform = meta.runtype == "ont" ? "--platform nanopore" : "--platform illumina"
    """
    tb-profiler \\
        profile \\
        $options.args \\
        $platform \\
        --csv \\
        --txt \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        --no_trim \\
        $input_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}
