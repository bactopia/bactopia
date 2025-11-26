nextflow.preview.types = true

process DEFENSEFINDER_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    db             : Path

    output:
    genes_tsv      = tuple(meta, files("*_defense_finder_genes.tsv"))
    hmmer_tsv      = tuple(meta, files("*_defense_finder_hmmer.tsv"))
    systems_tsv    = tuple(meta, files("*_defense_finder_systems.tsv"))
    proteins       = tuple(meta, files("*.prt"))
    proteins_index = tuple(meta, files("*.prt.idx"))
    macsydata_raw  = tuple(meta, file("${prefix}.macsydata.tar.gz", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
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
    """
    # Extract database
    # Use custom TMPDIR to prevent FileExistsError related to writing to same tmpdir (/tmp/tmp-macsy-cache/)
    tar -xf ${db}
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
        --workers ${task.cpus} \\
        --models-dir defense-finder/ \\
        ${fasta}

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
