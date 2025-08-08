process PHYLOFLASH  {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}/*")     , emit: results
    file "${prefix}/${prefix}.toalign.fasta"  , emit: aln, optional: true
    file "${prefix}/${prefix}.phyloFlash.json", emit: summary, optional: true
    path "versions.yml"                       , emit: versions
    path ".command.begin"                     , emit: begin
    path ".command.err"                       , emit: err
    path ".command.log"                       , emit: log
    path ".command.out"                       , emit: out
    path ".command.run"                       , emit: run
    path ".command.sh"                        , emit: sh
    path ".command.trace"                     , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
