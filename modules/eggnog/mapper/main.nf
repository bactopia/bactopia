nextflow.preview.types = true

process EGGNOG_MAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    db             : Path

    output:
    hits           = tuple(meta, file("*.emapper.hits"))
    seed_orthologs = tuple(meta, file("*.emapper.seed_orthologs"))
    annotations    = tuple(meta, file("*.emapper.annotations"))
    xlsx           = tuple(meta, file("*.emapper.annotations.xlsx", optional: true))
    orthologs      = tuple(meta, file("*.emapper.orthologs", optional: true))
    genepred       = tuple(meta, file("*.emapper.genepred.fasta", optional: true))
    gff            = tuple(meta, file("*.emapper.gff", optional: true))
    no_anno        = tuple(meta, file("*.emapper.no_annotations.fasta", optional: true))
    pfam           = tuple(meta, file("*.emapper.pfam", optional: true))
    logs           = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin       = tuple(meta, file(".command.begin"))
    nf_err         = tuple(meta, file(".command.err"))
    nf_log         = tuple(meta, file(".command.log"))
    nf_out         = tuple(meta, file(".command.out"))
    nf_run         = tuple(meta, file(".command.run"))
    nf_sh          = tuple(meta, file(".command.sh"))
    nf_trace       = tuple(meta, file(".command.trace"))
    versions       = tuple(meta, file("versions.yml"))

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
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        EGGNOG_DB=\$(find database/ -name "eggnog.db" | sed 's=eggnog.db==')
    else
        EGGNOG_DB=\$(find ${db}/ -name "eggnog.db" | sed 's=eggnog.db==')
    fi

    emapper.py \\
        ${task.ext.args} \\
        --cpu ${task.cpus} \\
        --data_dir \$EGGNOG_DB \\
        --output ${prefix} \\
        -i ${fasta}

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
