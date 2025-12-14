/**
 * Quality Assessment Tool for Genome Assemblies.
 *
 * Uses [QUAST](https://github.com/ablab/quast) to evaluate genome assemblies by computing various
 * metrics such as N50, gene counts, and assembly length.
 *
 * @status stable
 * @keywords quast, assembly, quality control, n50, metrics
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation quast
 *
 * @input tuple(meta, assembly, meta_file)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format (Path)
 * - `meta_file`: Meta file containing reference size information (Path)
 *
 * @output tsv          Transposed report in TSV format
 * @output supplemental Supplemental files including plots and HTML reports
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process QUAST {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly, meta_file) : Tuple<Map, Path, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, files("versions.yml"))

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

    est_ref_size=""
    # Use rev to get the last column easily, then re-reverse it
    ref_size=\$(tail -n 1 ${meta_file} | rev | cut -f 1 | rev)
    if [ "\${ref_size}" != "0" ]; then
        est_ref_size="--est-ref-size \${ref_size}"
    fi

    quast ${assembly_name} \${est_ref_size} \\
        -o supplemental \\
        --threads ${task.cpus} \\
        ${task.ext.args} \\
        --glimmer

    mv supplemental/quast.log ./
    mv supplemental/transposed_report.tsv ${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
