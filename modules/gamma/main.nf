/**
 * Identification, classification, and annotation of translated gene matches.
 *
 * Uses [GAMMA](https://github.com/rastanton/GAMMA) (Gene Allele Mutation Microbial Assessment)
 * to identify and annotate coding sequences in an assembly that match a specific gene database.
 * It is particularly useful for detecting specific targets like antimicrobial resistance genes
 * or virulence factors while accounting for potential mutations.
 *
 * @status stable
 * @keywords gene finding, annotation, homology, alignment, gamma, psl
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation gamma
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * The reference gene database in FASTA format
 *
 * @output gamma     The main GAMMA output file containing annotated gene matches
 * @output psl       The raw alignment details in PSL format
 * @output gff       Gene matches in GFF3 format
 * @output fasta     Extracted nucleotide sequences of the matched genes
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process GAMMA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>
    db             : Path

    output:
    gamma    = tuple(meta, files("*.gamma"))
    psl      = tuple(meta, files("*.psl"))
    gff      = tuple(meta, files("*.gff", optional: true))
    fasta    = tuple(meta, files("*.fasta", optional: true))
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
    def is_compressed = assembly.toList()[0].getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.toList()[0].getName().replace(".gz", "")
    def VERSION = '2.1'
    // Version information not provided by tool on CLI
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    GAMMA.py \\
        ${task.ext.args} \\
        ${assembly_name} \\
        ${db} \\
        ${prefix}

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: ${VERSION}
    END_VERSIONS
    """
}
