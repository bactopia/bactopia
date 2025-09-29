process CHECKM2_PREDICT {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: tsv
    tuple val(meta), path("supplemental/*"), emit: supplemental
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
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """    
    # Decompress fasta file if compressed
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Check if db is a directory - if so, find the diamond database
    if [ -d "$db" ]; then
        CHECKM2_DB=\$(find ${db}/ -name "*.dmnd")
    else
        CHECKM2_DB=$db
    fi

    checkm2 \\
        predict \\
        --output-directory supplemental \\
        --threads ${task.cpus} \\
        --database_path \$CHECKM2_DB \\
        $task.ext.args \\
        --input ${fasta}

    mv supplemental/checkm2.log ./
    mv supplemental/quality_report.tsv ./${prefix}.tsv

    # Cleanup
    gzip supplemental/protein_files/*.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
