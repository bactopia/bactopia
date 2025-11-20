nextflow.preview.types = true

process MLST {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    db             : Path

    output:
    tsv      = tuple(meta, file("*.tsv"))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
    versions = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    # Extract database
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        MLST_DB=\$(find database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    else
        MLST_DB=\$(find ${db}/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    fi

    mlst \\
        --threads ${task.cpus} \\
        --blastdb \$MLST_DB/blast/mlst.fa \\
        --datadir \$MLST_DB/pubmlst \\
        ${task.ext.args} \\
        ${fasta} \\
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
