process EGGNOG_MAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)
    path db

    output:
    tuple val(meta), path("*.emapper.hits")                , emit: hits
    tuple val(meta), path("*.emapper.seed_orthologs")      , emit: seed_orthologs
    tuple val(meta), path("*.emapper.annotations")         , emit: annotations
    tuple val(meta), path("*.emapper.annotations.xlsx")    , emit: xlsx     , optional: true
    tuple val(meta), path("*.emapper.orthologs")           , emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.genepred.fasta")      , emit: genepred , optional: true
    tuple val(meta), path("*.emapper.gff")                 , emit: gff      , optional: true
    tuple val(meta), path("*.emapper.no_annotations.fasta"), emit: no_anno  , optional: true
    tuple val(meta), path("*.emapper.pfam")                , emit: pfam     , optional: true
    tuple val(meta), path("*.{log,err}")                   , emit: logs     , optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        EGGNOG_DB=\$(find database/ -name "eggnog.db" | sed 's=eggnog.db==')
    else
        EGGNOG_DB=\$(find $db/ -name "eggnog.db" | sed 's=eggnog.db==')
    fi

    emapper.py \\
        ${task.ext.args} \\
        --cpu $task.cpus \\
        --data_dir \$EGGNOG_DB \\
        --output $prefix \\
        -i $fasta

    # Cleanup
    if [ "$is_tarball" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
