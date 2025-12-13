/**
 * Shigella and EIEC serotyping from assemblies.
 *
 * Uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to differentiate *Shigella* and
 * Enteroinvasive *E. coli* (EIEC) and predict their serotypes from genome assemblies. It utilizes
 * cluster-specific marker genes to distinguish these closely related pathovars.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, identification, cluster, virulence
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation shigeifinder
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv      ShigEiFinder results in TSV format
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SHIGEIFINDER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def VERSION = '1.3.2'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def is_compressed = assembly.toList()[0].getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.toList()[0].getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    shigeifinder \\
        ${task.ext.args} \\
        --output ${prefix}.tsv \\
        -t ${task.cpus} \\
        -i ${assembly_name}

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: ${VERSION}
    END_VERSIONS
    """
}
