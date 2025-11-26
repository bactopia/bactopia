nextflow.preview.types = true

process SHIGAPASS {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    flex_tsv = tuple(meta, file("${prefix}_Flex_summary.tsv", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    def SHIGAPASS_VERSION = "1.5.0"
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # ShigaPass does not accept compressed FASTA files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # Convert our genome path to a file with a path in it
    ls ${fasta_name} > ${fasta_name}_tmp.txt

    # Run ShigaPass
    ShigaPass.sh \\
        -l ${fasta_name}_tmp.txt \\
        ${task.ext.args} \\
        -p "\$(dirname \$(which ShigaPass.sh))/../share/shigapass-${SHIGAPASS_VERSION}/db" \\
        -t ${task.cpus} \\
        -o ${prefix}

    # Remove the temporary file from above
    rm ${fasta_name}_tmp.txt ${fasta_name}

    # Convert to tab delimited and move to the pwd
    sed 's/;/\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv

    # Convert to tab delimited and move to the pwd
    [ ! -f ${prefix}/ShigaPass_Flex_summary.csv ] || sed 's/;/\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """
}
