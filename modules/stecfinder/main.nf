/**
 * Find STEC gene markers in E. coli genomes.
 *
 * This process executes stecfinder to perform analysis
 *
 * @status stable
 * @keywords STEC, E. coli, virulence, typing
 * @tags complexity:moderate input-type:single output-type:single features:archive-output, compression, conditional-logic
 * @citation stecfinder
 *
 * @input tuple(meta, fasta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Assembly file in FASTA format
 * - `reads`: Reads file (if using reads mode)
 *
 * @output tsv      TSV file with STEC gene markers results
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process STECFINDER {
    tag "${prefix}"
    label 'process_low'
    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta, reads) : Tuple<Map, Path, List<Path>>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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
    def is_compressed = (fasta.getName().endsWith(".gz") ? true : false) && !task.ext.stecfinder_use_reads ? true : false
    def seq_name = is_compressed ? fasta.getName().replace(".gz", "") : reads.join(" ")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${seq_name}
    fi

    stecfinder \\
        -i ${seq_name} \\
        ${task.ext.args} \\
        -t ${task.cpus} > ${prefix}.tsv

    # Cleanup
    rm -rf ${seq_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stecfinder: \$(echo \$(stecfinder --version 2>&1) | sed 's/^.*STECFinder version: //;' )
    END_VERSIONS
    """
}
