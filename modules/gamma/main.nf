/**
 * Gene Allele Mutation Microbial Assessment.
 *
 * This process executes gamma to perform analysis
 *
 * @status stable
 * @keywords gamma, gene-calling
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic
 * @citation gamma
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: FASTA file
 *
 * @input db
 * Database in FASTA format
 *
 * @output gamma    GAMMA file with annotated gene matches
 * @output psl      PSL file with all gene matches found
 * @output gff      GFF file
 * @output fasta    multifasta file of the gene matches
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process GAMMA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    db             : Path

    output:
    gamma    = tuple(meta, files("*.gamma"))
    psl      = tuple(meta, files("*.psl"))
    gff      = tuple(meta, files("*.gff", optional: true))
    fasta    = tuple(meta, files("*.fasta", optional: true))
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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def VERSION = '2.1'
    // Version information not provided by tool on CLI
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    GAMMA.py \\
        ${task.ext.args} \\
        ${fasta_name} \\
        ${db} \\
        ${prefix}

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: ${VERSION}
    END_VERSIONS
    """
}
