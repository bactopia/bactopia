/**
 * Serovar prediction of Salmonella assemblies.
 *
 * Uses [SISTR](https://github.com/phac-nml/sistr_cmd) (Salmonella In Silico Typing Resource) to
 * predict serovars of *Salmonella* from draft genome assemblies using core genome Multi-Locus
 * Sequence Typing (cgMLST).
 *
 * @status stable
 * @keywords salmonella, serotype, cgmlst, typing, prediction
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation sistr
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv          SISTR prediction results in TSV format
 * @output allele_fasta Novel alleles in FASTA format
 * @output allele_json  Alleles in JSON format
 * @output cgmlst_csv   cgMLST profile in CSV format
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SISTR {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    allele_fasta = tuple(meta, file("${prefix}-allele.fasta.gz"))
    allele_json  = tuple(meta, file("${prefix}-allele.json.gz"))
    cgmlst_csv   = tuple(meta, file("${prefix}-cgmlst.csv"))
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

    sistr \\
        --qc \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        --alleles-output ${prefix}-allele.json \\
        --novel-alleles ${prefix}-allele.fasta \\
        --cgmlst-profiles ${prefix}-cgmlst.csv \\
        --output-prediction ${prefix} \\
        --output-format tab \\
        ${assembly_name}

    mv ${prefix}.tab ${prefix}.tsv
    gzip ${prefix}-allele.json
    gzip ${prefix}-allele.fasta

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sistr: \$(echo \$(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' )
    END_VERSIONS
    """
}
