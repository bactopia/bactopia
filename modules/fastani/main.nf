nextflow.preview.types = true

process FASTANI {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, query) : Tuple<Map, Set<Path>>
    reference      : Path

    stage:
    stageAs 'query-tmp/*', query

    output:
    tsv      = tuple(meta, file("*.tsv"))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
    versions = tuple(meta, file("versions.yml"))

    script:
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    reference_fasta = reference.getName().replace(".gz", "")
    reference_name = reference_fasta.replace(".fna", "")
    prefix = task.ext.prefix ?: "${reference_name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "fastani/logs/${prefix}"
    meta.process_name = task.ext.process_name
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${reference_fasta}
    fi

    mkdir query
    cp -L query-tmp/* query/
    find query/ -name "*.gz" | xargs gunzip
    find query/ -name "*" -type f > query-list.txt

    fastANI \\
        --ql query-list.txt \\
        -r ${reference_fasta} \\
        -o fastani-result.tmp

    echo "query<TAB>reference<TAB>ani<TAB>mapped_fragments<TAB>total_fragments" | sed 's/<TAB>/\t/g' > ${reference_name}.tsv
    sed 's=^query/==' fastani-result.tmp >> ${reference_name}.tsv

    # Cleanup
    rm -rf ${reference_fasta} query/ query-list.txt fastani-result.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
