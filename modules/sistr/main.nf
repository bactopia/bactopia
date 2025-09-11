process SISTR {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")            , emit: tsv
    tuple val(meta), path("*-allele.fasta.gz"), emit: allele_fasta
    tuple val(meta), path("*-allele.json.gz") , emit: allele_json
    tuple val(meta), path("*-cgmlst.csv")     , emit: cgmlst_csv
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
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    sistr \\
        --qc \\
        $args \\
        --threads $task.cpus \\
        --alleles-output ${prefix}-allele.json \\
        --novel-alleles ${prefix}-allele.fasta \\
        --cgmlst-profiles ${prefix}-cgmlst.csv \\
        --output-prediction ${prefix} \\
        --output-format tab \\
        $fasta_name

    mv ${prefix}.tab ${prefix}.tsv
    gzip ${prefix}-allele.json
    gzip ${prefix}-allele.fasta

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sistr: \$(echo \$(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' )
    END_VERSIONS
    """
}
