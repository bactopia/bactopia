process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(json, stageAs: 'results-tmp/*')

    output:
    tuple val(meta), path("tbprofiler.csv")         , emit: csv
    tuple val(meta), path("tbprofiler.variants.csv"), emit: variants_csv
    tuple val(meta), path("tbprofiler.variants.txt"), emit: variants_txt
    tuple val(meta), path("*.itol.*.txt")           , emit: itol, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    # Uncompress the JSON files
    mkdir results
    cp -L results-tmp/* results/
    find results/ -name "*.json.gz" | xargs gunzip

    tb-profiler \\
        collate \\
        ${task.ext.args} \\
        --db_dir database/ \\
        --format csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler collate --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
