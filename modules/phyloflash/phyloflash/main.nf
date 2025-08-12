process PHYLOFLASH  {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}/*")     , emit: results
    file "${prefix}/${prefix}.toalign.fasta"  , emit: aln, optional: true
    file "${prefix}/${prefix}.phyloFlash.json", emit: summary, optional: true
    tuple val(meta), path("*.{log,err}")      , emit: logs, optional: true
    tuple val(meta), path(".command.begin")   , emit: nf_begin
    tuple val(meta), path(".command.err")     , emit: nf_err
    tuple val(meta), path(".command.log")     , emit: nf_log
    tuple val(meta), path(".command.out")     , emit: nf_out
    tuple val(meta), path(".command.run")     , emit: nf_run
    tuple val(meta), path(".command.sh")      , emit: nf_sh
    tuple val(meta), path(".command.trace")   , emit: nf_trace
    tuple val(meta), path("versions.yml")     , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def read_opts = meta.single_end ? "-read1 ${reads[0]}" : "-read1 ${reads[0]} -read2 ${reads[1]}"
    """
    mkdir $prefix
    phyloFlash.pl \\
        $args \\
        $read_opts \\
        -lib $prefix \\
        -dbhome . \\
        -CPUs $task.cpus

    jsonify-phyloflash.py ${prefix}.phyloFlash > ${prefix}.phyloFlash.json
    mv ${prefix}.* $prefix


    if phyloflash-summary.py ${prefix}/ | grep -q -c "WARNING: Multiple SSUs were assembled by SPAdes"; then
        MULTI="1"
    fi

    if [ "${params.allow_multiple_16s}" == "true" ]; then
        MULTI="0"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloFlash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
