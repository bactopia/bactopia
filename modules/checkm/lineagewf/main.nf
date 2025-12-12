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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv           A tab-delimited report of quality metrics (Completeness, Contamination, Heterogeneity). See [CheckM's lineage_wf documentation](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) for details.
 * @output supplemental  Directory containing lineage files, marker gene stats, and storage logs
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process CHECKM_LINEAGEWF {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

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
