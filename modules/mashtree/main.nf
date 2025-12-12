/**
 * Rapid alignment-free phylogenomic tree construction.
 *
 * Uses [Mashtree](https://github.com/lskatz/mashtree) to create a phylogenetic tree
 * from genome sequences (FASTA, FASTQ, or GenBank) using MinHash distances. It computes
 * pairwise distances between all inputs and uses the Neighbor-Joining algorithm to
 * cluster genomes, effectively creating a "distance-based" tree without full alignment.
 *
 * @status stable
 * @keywords phylogeny, tree, mash, minhash, alignment-free, distance, clustering, neighbor-joining
 * @tags complexity:moderate input-type:multiple output-type:single features:conditional-logic
 * @citation mashtree
 *
 * @input tuple(meta, seqs)
 * - `meta`: Groovy Map containing sample information
 * - `seqs`: Assembled contigs in FASTA format
 *
 * @output tree      The final phylogenetic tree in Newick format (*.dnd)
 * @output matrix    The pairwise distance matrix used to build the tree (*.tsv)
 * @output sketches  Directory containing the individual Mash sketches (optional)
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MASHTREE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, seqs) : Tuple<Map, Set<Path>>

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
