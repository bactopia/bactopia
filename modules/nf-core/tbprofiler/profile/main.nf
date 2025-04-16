// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'tbprofiler')
options.btype = "tools"
conda_tools   = "bioconda::tb-profiler=6.6.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:6.6.3--pyhdfd78af_0' :
        'quay.io/biocontainers/tb-profiler:6.6.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bam/*.bam")        , emit: bam
    tuple val(meta), path("results/*.csv")    , emit: csv, optional: true
    tuple val(meta), path("results/*.json.gz"), emit: json
    tuple val(meta), path("results/*.txt")    , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")     , emit: vcf
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    def platform = meta.runtype == "ont" ? "--platform nanopore" : "--platform illumina"
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    tb-profiler \\
        profile \\
        $options.args \\
        $platform \\
        --csv \\
        --txt \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        --no_trim \\
        --db_dir database/ \\
        $input_reads

    # Cleanup
    gzip results/*.json

    tb-profiler profile --version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler profile --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
