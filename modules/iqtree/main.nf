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
 * @input record(meta, aln)
 * - `meta`: Groovy Record containing sample information
 * - `aln`: Multiple sequence alignment in FASTA, PHYLIP, or NEXUS format
 *
 * @output record(meta, aln, nwk, results, logs, nf_logs, versions)
 * - `aln`: The input alignment (passed through)
 * - `nwk`: The final maximum-likelihood phylogenetic tree (Newick format)
 *
 * @results iqtree (or iqtree-fast)
 * - `${prefix}.*`: IQ-TREE output files (model info, bootstrap trees, log, alninfo, etc.)
 */
nextflow.enable.types = true

// bactopia-lint: ignore M026
process IQTREE {
    tag "${prefix}"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        aln: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        aln: file("${aln}"),
        nwk: file(treefile),
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
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name == "iqtree-fast" ? "iqtree-fast" : task.ext.process_name
    args = process_name == "iqtree-fast" ? task.ext.fast_args : task.ext.args
    treefile = process_name == "iqtree-fast" ? "${process_name}/${prefix}.treefile" : "${prefix}.treefile"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "",
        logs_dir: "${process_name}/logs/",
        process_name: process_name
    )
    """
    iqtree \\
        ${args} \\
        -s ${aln} \\
        -nt ${task.cpus} \\
        -ntmax ${task.cpus} \\
        -pre ${prefix}

    # Only gzip files if they exist
    if [[ -f "${prefix}.alninfo" ]]; then
        gzip ${prefix}.alninfo
    fi

    mkdir ${process_name}/
    find . -maxdepth 1 -name "${prefix}*" -not -name "${aln}" -not -type d -exec mv {} ${process_name}/ \\;

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
