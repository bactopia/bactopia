process MYKROBE_PREDICT {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(seqs)
    val species

    output:
    tuple val(meta), path("${prefix}.csv") , emit: csv
    tuple val(meta), path("${prefix}.json"), emit: json
    path "versions.yml"                    , emit: versions
    path ".command.begin"                  , emit: begin
    path ".command.err"                    , emit: err
    path ".command.log"                    , emit: log
    path ".command.out"                    , emit: out
    path ".command.run"                    , emit: run
    path ".command.sh"                     , emit: sh
    path ".command.trace"                  , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    is_ont = meta.runtype == "ont" ? "--ont" : ""
    """
    mykrobe \\
        predict \\
        $args $is_ont \\
        --species $species \\
        --threads $task.cpus \\
        --sample $prefix \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq $seqs

    # Cleanup
    rm -rf mykrobe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
