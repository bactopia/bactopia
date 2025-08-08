process DEFENSEFINDER_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)
    each path(db)

    output:
    tuple val(meta), path("*_defense_finder_genes.tsv")  , emit: genes_tsv
    tuple val(meta), path("*_defense_finder_hmmer.tsv")  , emit: hmmer_tsv
    tuple val(meta), path("*_defense_finder_systems.tsv"), emit: systems_tsv
    tuple val(meta), path("*.prt")                       , emit: proteins
    tuple val(meta), path("*.prt.idx")                   , emit: proteins_index
    tuple val(meta), path("${prefix}.macsydata.tar.gz")  , emit: macsydata_raw, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
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
        $args \\
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
