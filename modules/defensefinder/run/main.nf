process DEFENSEFINDER_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)
    each path(db)

    output:
    tuple val(meta), path("*_defense_finder_genes.tsv")  , emit: genes_tsv
    tuple val(meta), path("*_defense_finder_hmmer.tsv")  , emit: hmmer_tsv
    tuple val(meta), path("*_defense_finder_systems.tsv"), emit: systems_tsv
    tuple val(meta), path("*.prt")                       , emit: proteins
    tuple val(meta), path("*.prt.idx")                   , emit: proteins_index
    tuple val(meta), path("${prefix}.macsydata.tar.gz")  , emit: macsydata_raw, optional: true
    tuple val(meta), path("*.{log,err}")                 , emit: logs, optional: true
    tuple val(meta), path(".command.begin")              , emit: nf_begin
    tuple val(meta), path(".command.err")                , emit: nf_err
    tuple val(meta), path(".command.log")                , emit: nf_log
    tuple val(meta), path(".command.out")                , emit: nf_out
    tuple val(meta), path(".command.run")                , emit: nf_run
    tuple val(meta), path(".command.sh")                 , emit: nf_sh
    tuple val(meta), path(".command.trace")              , emit: nf_trace
    tuple val(meta), path("versions.yml")                , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    # Extract database
    # Use custom TMPDIR to prevent FileExistsError related to writing to same tmpdir (/tmp/tmp-macsy-cache/)
    tar -xf $db
    mkdir -p df-tmp/df
    TMPDIR=df-tmp/df HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/defense-finder-models-v${task.ext.df_models_version}.tar.gz

    mkdir -p df-tmp/cf
    TMPDIR=df-tmp/cf HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/CasFinder-${task.ext.casfinder_version}.tar.gz

    TMPDIR=df-tmp/ HOME=df-tmp/ defense-finder \\
        run \\
        ${task.ext.args} \\
        --workers $task.cpus \\
        --models-dir defense-finder/ \\
        $fasta

    if [ "${task.ext.df_preserveraw}" == "true" ]; then
        tar -czf ${prefix}.macsydata.tar.gz defense-finder-tmp/
        rm -rf defense-finder-tmp/ 
    fi

    # Cleanup intermediate files and unused outputs
    rm -rf models/ defense-finder/ df-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${task.ext.df_version}
        defense-finder-models: ${task.ext.df_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
