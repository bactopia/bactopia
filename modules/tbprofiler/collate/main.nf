process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda_env}"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:6.6.3--pyhdfd78af_0' :
        'quay.io/biocontainers/tb-profiler:6.6.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(json, stageAs: 'results-tmp/*')

    output:
    tuple val(meta), path("tbprofiler.csv")   , emit: csv
    tuple val(meta), path("tbprofiler.variants.csv"), emit: variants_csv
    tuple val(meta), path("tbprofiler.variants.txt"), emit: variants_txt
    tuple val(meta), path("*.itol.*.txt")            , emit: itol, optional: true
    path "*.{log,err}"                              , emit: logs, optional: true
    path ".command.begin"   , emit: begin
    path ".command.err"     , emit: err
    path ".command.log"     , emit: log
    path ".command.out"     , emit: out
    path ".command.run"     , emit: run
    path ".command.sh"      , emit: sh
    path ".command.trace"   , emit: trace
    path "versions.yml"                             , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
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
