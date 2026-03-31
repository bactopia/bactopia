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
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, nwk, tsv, sketches, results, logs, nf_logs, versions)
 * - `nwk`: The final phylogenetic tree in Newick format (*.dnd)
 * - `tsv`: The pairwise distance matrix used to build the tree (*.tsv)
 * - `sketches`: Directory containing the individual Mash sketches (optional)
 */
nextflow.preview.types = true

process MASHTREE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Set<Path>): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        nwk: file("${prefix}.dnd"),
        tsv: file("${prefix}.tsv"),
        sketches: files("sketches/*", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.dnd"),
            files("${prefix}.tsv"),
            files("sketches/*", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

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
        ${fna.join(' ')}

    # Cleanup
    rm -rf mashtree-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
