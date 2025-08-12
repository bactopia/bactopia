process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bam/*.bam")        , emit: bam
    tuple val(meta), path("results/*.csv")    , emit: csv, optional: true
    tuple val(meta), path("results/*.json.gz"), emit: json
    tuple val(meta), path("results/*.txt")    , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.{log,err}")      , emit: logs, optional: true
    tuple val(meta), path(".command.begin")   , emit: nf_begin
    tuple val(meta), path(".command.err")     , emit: nf_err
    tuple val(meta), path(".command.log")     , emit: nf_log
    tuple val(meta), path(".command.out")     , emit: nf_out
    tuple val(meta), path(".command.run")     , emit: nf_run
    tuple val(meta), path(".command.sh")      , emit: nf_sh
    tuple val(meta), path(".command.trace")   , emit: nf_trace
    tuple val(meta), path("versions.yml")     , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
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
