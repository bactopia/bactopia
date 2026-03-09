/**
 * Assess genome quality using lineage-specific marker sets.
 *
 * Uses [CheckM](https://github.com/Ecogenomics/CheckM) to estimate the completeness and
 * contamination of genome assemblies. It places the genome into a reference tree to select
 * an appropriate set of single-copy marker genes, then calculates quality metrics based on
 * the recovery of these markers.
 *
 * @status stable
 * @keywords quality control, completeness, contamination, marker genes, lineage, bacteria, archaea
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation checkm
 *
 * @note Database Required
 * Requires the CheckM reference database (~275GB uncompressed) to be configured via the
 * `CHECKM_DATA_PATH` environment variable or pre-installed in the container.
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 */
nextflow.preview.types = true

process CHECKM_LINEAGEWF {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        tsv: file("${prefix}.tsv"),
        results: files("${prefix}.tsv") + files("supplemental/*"),
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
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

    checkm \\
        lineage_wf ./ supplemental/ \\
        --tab_table \\
        --threads ${task.cpus} \\
        --pplacer_threads ${task.cpus} \\
        --alignment_file supplemental/${prefix}-genes.aln \\
        --file supplemental/${prefix}-results.txt \\
        ${task.ext.args}

    find ./supplemental/ -name "*.faa" -or -name "*hmmer.analyze.txt" -or -name "*.fasta" | xargs gzip
    mv supplemental/checkm.log ./
    mv supplemental/${prefix}-results.txt ./${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(echo \$(checkm -h 2>&1) | sed 's/.*CheckM v//;s/ .*\$//')
    END_VERSIONS
    """
}
