/**
 * Predict Penicillin Binding Protein (PBP) type of *Streptococcus pneumoniae* assemblies.
 *
 * Uses [PBPtyper](https://github.com/rpetit3/pbptyper) to detect variations in the three
 * key PBP genes (*pbp1a*, *pbp2b*, and *pbp2x*) in *S. pneumoniae*. Typing these genes is
 * essential for predicting reduced susceptibility or full resistance to penicillin and other
 * $\beta$-lactam antibiotics.
 *
 * @status stable
 * @keywords bacteria, streptococcus pneumoniae, penicillin, amr, resistance, pbp, typing
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation pbptyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, blast, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary file with the predicted PBP type for each gene
 * - `blast`: A tab-delimited file of the raw TBLASTN hits used for gene identification
 */
nextflow.preview.types = true

process PBPTYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        fna: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        blast: file("${prefix}.tblastn.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.tblastn.tsv")
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    pbptyper \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${fna}

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/.*pbptyper, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pbptyper --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
