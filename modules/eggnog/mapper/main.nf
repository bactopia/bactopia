process EGGNOG_MAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)
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
    path "*.{log,err}"                                     , emit: logs     , optional: true
    path ".command.{begin,err,log,out,run,sh,trace}"       , emit: nf_logs
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        EGGNOG_DB=\$(find database/ -name "eggnog.db" | sed 's=eggnog.db==')
    else
        EGGNOG_DB=\$(find $db/ -name "eggnog.db" | sed 's=eggnog.db==')
    fi

    emapper.py \\
        $args \\
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
