process MLST {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions
    path ".command.begin"         , emit: begin
    path ".command.err"           , emit: err
    path ".command.log"           , emit: log
    path ".command.out"           , emit: out
    path ".command.run"           , emit: run
    path ".command.sh"            , emit: sh
    path ".command.trace"         , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    # Extract database
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        MLST_DB=\$(find database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    else
        MLST_DB=\$(find $db/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    fi

    mlst \\
        --threads $task.cpus \\
        --blastdb \$MLST_DB/blast/mlst.fa \\
        --datadir \$MLST_DB/pubmlst \\
        $args \\
        $fasta \\
        > ${prefix}.tsv

    if [[ -f "\$MLST_DB/DB_VERSION" ]]; then
        DB_VERSION=\$(cat \$MLST_DB/DB_VERSION)
    else
        DB_VERSION="custom database"
    fi

    # Cleanup
    rm -rf database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/^.*mlst //' )
        database: \$DB_VERSION
    END_VERSIONS
    """
}
