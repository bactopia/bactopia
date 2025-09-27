process GTDBTK_CLASSIFYWF {
    tag "${prefix}"
    label 'process_high'
    label 'process_high_memory'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fna, stageAs: 'fna-tmp/*')
    path db, stageAs: 'gtdb/*'

    output:
    tuple val(meta), path("supplemental/*")         , emit: results
    tuple val(meta), path("${prefix}.*.summary.tsv"), emit: tsv
    tuple val(meta), path("*.{log,err}")            , emit: logs, optional: true
    tuple val(meta), path(".command.begin")         , emit: nf_begin
    tuple val(meta), path(".command.err")           , emit: nf_err
    tuple val(meta), path(".command.log")           , emit: nf_log
    tuple val(meta), path(".command.out")           , emit: nf_out
    tuple val(meta), path(".command.run")           , emit: nf_run
    tuple val(meta), path(".command.sh")            , emit: nf_sh
    tuple val(meta), path(".command.trace")         , emit: nf_trace
    tuple val(meta), path("versions.yml")           , emit: versions

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
        export GTDBTK_DATA_PATH="\$(realpath \$(find database/ -path "*metadata*" -name "metadata.txt" | sed 's=/metadata/metadata.txt=='))"
    else
        export GTDBTK_DATA_PATH="\$(readlink $db)"
    fi
    mkdir fna
    cp -L fna-tmp/* fna/
    find fna/ -name "*.fna.gz" | xargs gunzip

    gtdbtk classify_wf \\
        ${task.ext.args} \\
        --cpus $task.cpus \\
        --pplacer_cpus $task.cpus \\
        --genome_dir ./fna \\
        --out_dir supplemental \\
        --skip_ani_screen \\
        --prefix ${prefix}
    mv supplemental/*.log ./
    mv supplemental/*.summary.tsv ./

    # Cleanup
    if [ "$is_tarball" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi
    if [ "${task.ext.gtdb_keep_msa}" == "false" ]; then
        # Delete MSA of submitted and reference genomes.
        rm -rf supplemental/align/*.msa.fasta.gz
    fi
    rm -rf fna/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
