/**
 * Efficient phylogenomic inference using Maximum Likelihood.
 *
 * Uses [IQ-TREE](http://www.iqtree.org/) to construct a maximum-likelihood phylogenetic tree
 * from a multiple sequence alignment. It automatically determines the best-fit substitution model
 * (via ModelFinder) and assesses branch support using the Ultrafast Bootstrap approximation.
 *
 * @status stable
 * @keywords phylogeny, tree, maximum likelihood, bootstrap, model selection, iqtree
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic
 * @citation iqtree
 *
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA, PHYLIP, or NEXUS format
 *
 * @output supplemental  Directory containing the detailed report (*.iqtree), distance matrix, and model parameters
 * @output phylogeny     The final maximum-likelihood phylogenetic tree (Newick format)
 * @output alignment     The input alignment (passed through)
 * @output aln_tree      A convenience tuple containing both the alignment and the tree
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process IQTREE {
    tag "${prefix}"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, msa) : Tuple<Map, Path>

    output:
    phylogeny    = tuple(meta, file(treefile))
    aln_tree     = tuple(meta, msa, file(treefile))
    supplemental = tuple(meta, files("${process_name}/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, files("versions.yml"))

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

    mkdir temp
    mv ${prefix}* temp/
    mv temp/ ${process_name}/

    if [ "${process_name}" != "iqtree-fast" ]; then
        mv ${process_name}/${prefix}.treefile ./
        mv ${process_name}/${msa} ./
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
