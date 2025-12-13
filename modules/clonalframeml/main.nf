/**
 * Inference of recombination in bacterial genomes.
 *
 * Uses [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) to detect recombination
 * events in bacterial genomes. It corrects the phylogenetic tree for recombination and produces
 * a "masked" alignment where recombinant regions are removed, allowing for more accurate
 * phylogenetic inference.
 *
 * @status stable
 * @keywords bacteria, recombination, phylogeny, alignment, msa, evolution
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic,compression
 * @citation clonalframeml
 *
 * @input tuple(meta, msa, newick)
 * - `meta`: Groovy Map containing sample information
 * - `msa`: Multiple sequence alignment in FASTA format
 * - `newick`: Initial phylogenetic tree in Newick format
 *
 * @output emsim       Uncertainty estimation results (if requested)
 * @output em          Final parameter estimates from the EM algorithm
 * @output status      Tab-delimited list of predicted recombination events (importations)
 * @output newick      The input tree with internal nodes labelled
 * @output fasta       Reconstructed ancestral sequences (*.fasta.gz)
 * @output pos_ref     Position cross-reference table (*.txt.gz)
 * @output masked_aln  The input alignment with recombinant regions masked (*.aln.gz)
 * @output logs        Optional software execution logs containing warnings/errors
 * @output nf_logs     Nextflow execution scripts and logs for debugging
 * @output versions    A YAML formatted file with software versions
 */
nextflow.preview.types = true

process CLONALFRAMEML {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, msa, newick) : Tuple<Map, Path, Path>

    output:
    emsim      = tuple(meta, files("${task.ext.process_name}/${prefix}.emsim.txt", optional: true))
    em         = tuple(meta, files("${task.ext.process_name}/${prefix}.em.txt"))
    status     = tuple(meta, files("${task.ext.process_name}/${prefix}.importation_status.txt"))
    newick     = tuple(meta, files("${task.ext.process_name}/${prefix}.labelled_tree.newick"))
    fasta      = tuple(meta, files("${task.ext.process_name}/${prefix}.ML_sequence.fasta.gz"))
    pos_ref    = tuple(meta, files("${task.ext.process_name}/${prefix}.position_cross_reference.txt.gz"))
    masked_aln = tuple(meta, files("${prefix}.masked.aln.gz"))
    logs       = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs    = tuple(meta, files(".command.*"))
    versions   = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.ext.process_name}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${task.ext.process_name}/logs/"
    meta.process_name = task.ext.process_name

    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${msa} > ${msa_name}
    fi

    ClonalFrameML \\
        ${newick} \\
        ${msa_name} \\
        ${prefix} \\
        ${task.ext.args}

    maskrc-svg.py ${prefix} --aln ${msa_name} --symbol '-' --out ${prefix}.masked.aln
    gzip ${prefix}.masked.aln

    # Cleanup
    rm ${msa_name}
    gzip *.ML_sequence.fasta
    gzip *.position_cross_reference.txt

    # Organize 
    mkdir ${task.ext.process_name}/
    mv *.emsim.txt ${task.ext.process_name}/
    mv *.em.txt ${task.ext.process_name}/
    mv *.importation_status.txt ${task.ext.process_name}/
    mv *.labelled_tree.newick ${task.ext.process_name}/
    mv *.ML_sequence.fasta.gz ${task.ext.process_name}/
    mv *.position_cross_reference.txt.gz ${task.ext.process_name}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clonalframeml: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
        maskrc-svg: \$( echo \$(maskrc-svg.py --version 2>&1) | sed 's/^.*maskrc-svg.py //' )
    END_VERSIONS
    """
}
