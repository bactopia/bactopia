process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda_env}"
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
    path ".command.begin"   , emit: begin
    path ".command.err"     , emit: err
    path ".command.log"     , emit: log
    path ".command.out"     , emit: out
    path ".command.run"     , emit: run
    path ".command.sh"      , emit: sh
    path ".command.trace"   , emit: trace
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    def platform = meta.runtype == "ont" ? "--platform nanopore" : "--platform illumina"
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    tb-profiler \\
        profile \\
        ${task.ext.args} \\
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
