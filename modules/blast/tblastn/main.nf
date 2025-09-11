process BLAST_TBLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(blastdb)
    path query

    output:
    tuple val(meta), path('*.tblastn.tsv') , emit: tsv
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
    // genes -> ffn, contigs -> fna
    def db_type = task.ext.use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf $blastdb
    
    ${which_cat} $query | \\
    tblastn \\
        -num_threads $task.cpus \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        $task.ext.args \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "$outcols" | sed 's/<TAB>/\t/g' > ${prefix}.tblastn.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastn.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
