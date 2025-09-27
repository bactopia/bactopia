process HICAP {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)
    path database_dir
    path model_fp

    output:
    tuple val(meta), path("*.gbk"), emit: gbk, optional: true
    tuple val(meta), path("*.svg"), emit: svg, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv
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
    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    hicap \\
        --query_fp $fasta_name \\
        $database_args \\
        $model_args \\
        ${task.ext.args} \\
        --threads $task.cpus \\
        --debug \\
        -o ./

    if [ ! -f ${prefix}.tsv ]; then
        echo "isolate<TAB>predicted_serotype<TAB>attributes<TAB>genes_identified<TAB>locus_location<TAB>region_I_genes<TAB>region_II_genes<TAB>region_III_genes<TAB>IS1016_hits" | sed 's/<TAB>/\t/g' > ${prefix}.tsv
        echo "${prefix}<TAB>cap_not_found<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-" | sed 's/<TAB>/\t/g' >> ${prefix}.tsv
    else
        sed -i 's/#isolate/isolate/' ${prefix}.tsv
    fi

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
