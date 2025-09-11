process ISMAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    path(reference)
    path(query)

    output:
    tuple val(meta), path("supplemental/*"), emit: results
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def ref_compressed = reference.getName().endsWith(".gz") ? true : false
    def reference_name = reference.getName().replace(".gz", "")
    def query_compressed = query.getName().endsWith(".gz") ? true : false
    def query_name = query.getName().replace(".gz", "")
    """
    if [ "$ref_compressed" == "true" ]; then
        gzip -c -d $reference > $reference_name
    fi
    if [ "$query_compressed" == "true" ]; then
        gzip -c -d $query > $query_name
    fi
    
    ismap \\
        $task.ext.args \\
        --t $task.cpus \\
        --output_dir $prefix \\
        --queries $query_name \\
        --log ${prefix} \\
        --reference $reference_name \\
        --reads $reads

    # Reorganize output files
    mkdir supplemental
    mv $prefix/*/* supplemental/

    # Cleanup and compress FASTQ and BED files
    rm -rf ${reference_name} ${query_name} ${prefix}/
    find supplemental/ -name "*.fastq" | xargs -I {} gzip {}
    find supplemental/ -name "*.bed" | xargs -I {} gzip {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ismapper: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
