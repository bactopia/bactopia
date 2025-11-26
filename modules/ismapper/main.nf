nextflow.preview.types = true

process ISMAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, Path>
    reference      : Path
    query          : Path

    output:
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    def query_name = query.getName().replace(".gz", "")
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${query_name}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${query_name}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def ref_compressed = reference.getName().endsWith(".gz") ? true : false
    def reference_name = reference.getName().replace(".gz", "")
    def query_compressed = query.getName().endsWith(".gz") ? true : false
    """
    if [ "${ref_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${reference_name}
    fi
    if [ "${query_compressed}" == "true" ]; then
        gzip -c -d ${query} > ${query_name}
    fi
    
    ismap \\
        ${task.ext.args} \\
        --t ${task.cpus} \\
        --output_dir ${prefix} \\
        --queries ${query_name} \\
        --log ${prefix} \\
        --reference ${reference_name} \\
        --reads ${reads}

    # Reorganize output files
    mkdir supplemental
    mv ${prefix}/*/* supplemental/

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
