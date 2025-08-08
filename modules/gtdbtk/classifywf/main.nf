process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fna, stageAs: 'fna-tmp/*')
    path db, stageAs: 'gtdb/*'

    output:
    path "results/*"                                        , emit: results
    tuple val(meta), path("results/${prefix}.*.summary.tsv"), emit: tsv
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml", emit: versions

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
        export GTDBTK_DATA_PATH="\$(realpath \$(find database/ -path "*metadata*" -name "metadata.txt" | sed 's=/metadata/metadata.txt=='))"
    else
        export GTDBTK_DATA_PATH="\$(readlink $db)"
    fi
    mkdir fna
    cp -L fna-tmp/* fna/
    find fna/ -name "*.fna.gz" | xargs gunzip

    gtdbtk classify_wf \\
        $args \\
        --cpus $task.cpus \\
        --pplacer_cpus $task.cpus \\
        --genome_dir ./fna \\
        --out_dir results \\
        --skip_ani_screen \\
        --prefix ${prefix}
    mv results/*.log ./

    # Cleanup
    if [ "$is_tarball" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi
    if [ "${task.ext.gtdb_keep_msa}" == "false" ]; then
        # Delete MSA of submitted and reference genomes.
        rm -rf results/align/*.msa.fasta.gz
    fi
    rm -rf fna/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
