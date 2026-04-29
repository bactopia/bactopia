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
 * @citation clonalframeml, maskrc_svg
 *
 * @input record(meta, aln, nwk)
 * - `meta`: Groovy Record containing sample information
 * - `aln`: Multiple sequence alignment in FASTA format
 * - `nwk`: Initial phylogenetic tree in Newick format
 *
 * @output record(meta, emsim?, em, status, nwk, fasta, pos_ref, masked_aln, results, logs, nf_logs, versions)
 * - `emsim?`: Uncertainty estimation results (if requested)
 * - `em`: Final parameter estimates from the EM algorithm
 * - `status`: Tab-delimited list of predicted recombination events (importations)
 * - `nwk`: The input tree with internal nodes labelled
 * - `fasta`: Reconstructed ancestral sequences (*.fasta.gz)
 * - `pos_ref`: Position cross-reference table (*.txt.gz)
 * - `masked_aln`: The input alignment with recombinant regions masked (*.aln.gz)
 */
nextflow.enable.types = true

process CLONALFRAMEML {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        aln: Path,
        nwk: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        emsim: file("${task.ext.process_name}/${prefix}.emsim.txt", optional: true),
        em: file("${task.ext.process_name}/${prefix}.em.txt"),
        status: file("${task.ext.process_name}/${prefix}.importation_status.txt"),
        nwk: file("${task.ext.process_name}/${prefix}.labelled_tree.newick"),
        fasta: file("${task.ext.process_name}/${prefix}.ML_sequence.fasta.gz"),
        pos_ref: file("${task.ext.process_name}/${prefix}.position_cross_reference.txt.gz"),
        masked_aln: file("${prefix}.masked.aln.gz"),
        // Generic fields (used for publishing)
        results: [
            files("${task.ext.process_name}/${prefix}.emsim.txt", optional: true),
            files("${task.ext.process_name}/${prefix}.em.txt"),
            files("${task.ext.process_name}/${prefix}.importation_status.txt"),
            files("${task.ext.process_name}/${prefix}.labelled_tree.newick"),
            files("${task.ext.process_name}/${prefix}.ML_sequence.fasta.gz"),
            files("${task.ext.process_name}/${prefix}.position_cross_reference.txt.gz"),
            files("${prefix}.masked.aln.gz")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.ext.process_name}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "",
        logs_dir: "${task.ext.process_name}/logs/",
        process_name: task.ext.process_name
    )

    def is_compressed = aln.getName().endsWith(".gz") ? true : false
    def aln_name = aln.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${aln} > ${aln_name}
    fi

    ClonalFrameML \\
        ${nwk} \\
        ${aln_name} \\
        ${prefix} \\
        ${task.ext.args}

    maskrc-svg.py ${prefix} --aln ${aln_name} --symbol '-' --out ${prefix}.masked.aln
    gzip ${prefix}.masked.aln

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm ${aln_name}
    fi
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
