/**
 * Quickly create a tree using Mash distances.
 *
 * This process executes mashtree to perform analysis
 *
 * @status stable
 * @keywords tree, mash, fasta, fastq
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation mashtree
 *
 * @input tuple(meta, seqs)
 * - `meta`: Groovy Map containing sample information
 * - `seqs`: FASTA, FASTQ, GenBank, or Mash sketch files
 *
 * @output tree     A Newick formatted tree file
 * @output matrix   A TSV matrix of pair-wise Mash distances
 * @output sketches Sketches
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process MASHTREE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, seqs) : Tuple<Map, Path>

    output:
    tree     = tuple(meta, file("${prefix}.dnd"))
    matrix   = tuple(meta, file("${prefix}.tsv"))
    sketches = tuple(meta, files("sketches/*", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "logs/"
    meta.process_name = task.ext.process_name
    """
    mkdir mashtree-tmp

    mashtree \\
        ${task.ext.args} \\
        --numcpus ${task.cpus} \\
        --outmatrix ${prefix}.tsv \\
        --outtree ${prefix}.dnd \\
        --tempdir mashtree-tmp/ \\
        ${seqs}

    # Clean up
    rm -rf mashtree-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
