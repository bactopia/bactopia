/**
 * Serotyping and Multi-Antigen Sequence Typing (MAST) of *Neisseria gonorrhoeae*.
 *
 * Uses [NG-MASTER](https://github.com/phac-nml/NG-MASTER) to identify the alleles of the
 * *porB* and *tbpB* genes in *N. gonorrhoeae* assemblies, which is the basis for the
 * internationally recognized MAST genotyping scheme.
 *
 * @status stable
 * @keywords bacteria, neisseria gonorrhoeae, serotype, typing, mast, porb, tbpb
 * @tags complexity:simple input-type:single output-type:single features:compression,conditional-logic
 * @citation ngmaster
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, results, logs, nf_logs, versions)
 * - `meta`: Groovy Map containing sample information and output paths
 * - `results`: List of result files including the TSV summary
 * - `logs`: Optional software execution logs containing warnings/errors
 * - `nf_logs`: Nextflow execution scripts and logs for debugging
 * - `versions`: A YAML formatted file with software versions
 */
nextflow.preview.types = true

process NGMASTER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta:     meta,
        // Named field (upstream consumers access this)
        tsv:      file("${prefix}.tsv"),
        // Generic fields (same convention across every module)
        results:  [file("${prefix}.tsv")],
        logs:     files("*.{log,err}", optional: true),
        nf_logs:  files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    ngmaster \\
        ${task.ext.args} \\
        ${assembly_name} \\
        > ${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngmaster: \$( echo \$(ngmaster --version 2>&1) | sed 's/^.*ngmaster //' )
    END_VERSIONS
    """
}
