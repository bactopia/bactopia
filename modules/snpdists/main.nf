/**
 * Create a SNP distance matrix from a multiple sequence alignment.
 *
 * Uses [snp-dists](https://github.com/tseemann/snp-dists) to read a FASTA alignment and
 * compute a pairwise SNP distance matrix between all sequences.
 *
 * @status stable
 * @keywords snp, distance, matrix, alignment, phylogeny
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation snpdists
 *
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output tsv      Pairwise SNP distance matrix in TSV format
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SNPDISTS {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, alignment) : Tuple<Map, Set<Path>>

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name ?: task.ext.process_name

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.ext.process_name}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${process_name}/logs/"
    meta.process_name = process_name
    """
    snp-dists \\
        ${task.ext.args} \\
        ${alignment} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-dists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
