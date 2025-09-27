process FASTANI {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(query, stageAs: 'query-tmp/*')
    each path(reference)

    output:
    tuple val(meta), path("*.tsv")         , emit: tsv
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
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    reference_fasta = reference.getName().replace(".gz", "")
    reference_name = reference_fasta.replace(".fna", "")
    prefix = task.ext.prefix ?: "${reference_name}"
    
    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${task.ext.rundir}/"
    meta.logs_dir = "${task.ext.rundir}/fastani/logs/${prefix}"
    meta.process_name = task.ext.process_name
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $reference_fasta
    fi

    mkdir query
    cp -L query-tmp/* query/
    find query/ -name "*.gz" | xargs gunzip
    find query/ -name "*" -type f > query-list.txt

    fastANI \\
        --ql query-list.txt \\
        -r $reference_fasta \\
        -o fastani-result.tmp

    echo "query<TAB>reference<TAB>ani<TAB>mapped_fragments<TAB>total_fragments" | sed 's/<TAB>/\t/g' > ${reference_name}.tsv
    sed 's=^query/==' fastani-result.tmp >> ${reference_name}.tsv

    # Cleanup
    rm -rf $reference_fasta query/ query-list.txt fastani-result.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
