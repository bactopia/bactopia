/**
 * Efficient phylogenomic inference using Maximum Likelihood.
 *
 * Uses [IQ-TREE](https://iqtree.github.io/) to construct a maximum-likelihood phylogenetic tree
 * from a multiple sequence alignment. It automatically determines the best-fit substitution model
 * (via ModelFinder) and assesses branch support using the Ultrafast Bootstrap approximation.
 *
 * @status stable
 * @keywords phylogeny, tree, maximum likelihood, bootstrap, model selection, iqtree
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic
 * @citation iqtree
 *
 * @input record(meta, msa)
 * - `meta`: Groovy Map containing sample information
 * - `msa`: Multiple sequence alignment in FASTA, PHYLIP, or NEXUS format
 *
 * @output record(meta, msa, phylogeny, results, logs, nf_logs, versions)
 * - `msa`: The input alignment (passed through)
 * - `phylogeny`: The final maximum-likelihood phylogenetic tree (Newick format)
 *
 * @results iqtree (or iqtree-fast)
 * - `${prefix}.*`: IQ-TREE output files (model info, bootstrap trees, log, alninfo, etc.)
 */
nextflow.preview.types = true

// bactopia-lint: ignore M026
process IQTREE {
    tag "${prefix}"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, msa: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        msa: file(msa),
        phylogeny: file(treefile),
        // Generic fields (used for publishing)
        results: [
            files(treefile),
            files("${process_name}/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name == "iqtree-fast" ? "iqtree-fast" : task.ext.process_name
    args = process_name == "iqtree-fast" ? task.ext.fast_args : task.ext.args
    treefile = process_name == "iqtree-fast" ? "${process_name}/${prefix}.treefile" : "${prefix}.treefile"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${process_name}/logs/"
    meta.process_name = process_name
    """
    iqtree \\
        ${args} \\
        -s ${msa} \\
        -nt ${task.cpus} \\
        -ntmax ${task.cpus} \\
        -pre ${prefix}

    # Only gzip files if they exist
    if [[ -f "${prefix}.alninfo" ]]; then
        gzip ${prefix}.alninfo
    fi

    mkdir ${process_name}/
    find . -maxdepth 1 -name "${prefix}*" -not -name "${msa}" -not -type d -exec mv {} ${process_name}/ \\;

    if [ "${process_name}" == "iqtree" ]; then
        # We don't want the fast-tree to be on the same level as the main tree in the outputs
        mv ${process_name}/${prefix}.treefile ./
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE version //;s/ .*//')
    END_VERSIONS
    """
}
