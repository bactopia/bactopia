/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * Screens assemblies for antimicrobial resistance and virulence genes using
 * [Abricate](https://github.com/tseemann/abricate). It bundles several databases
 * including NCBI, CARD, ResFinder, PlasmidFinder, ARG-ANNOT, and VFDB.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, antimicrobial resistance, virulence, plasmid, mobile genetic elements
 * @tags complexity:simple input-type:single output-type:single
 * @citation abricate
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output report   A tab-delimited report of hits, for full details please see [Abricate - Output](https://github.com/tseemann/abricate#output)
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process ABRICATE_RUN {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker }"

    input:
    (_meta, assembly): Tuple<Map, Set<Path>>

    output:
    report   = tuple(meta, files("${prefix}.tsv"))
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
    """
    abricate \\
        $assembly \\
        $task.ext.args \\
        --threads $task.cpus > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
