/**
 * In silico typing of the H. influenzae capsule locus.
 *
 * This process executes hicap to perform analysis
 *
 * @status stable
 * @keywords haemophilus, influenzae, capsule, typing, serotype
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent
 * @citation hicap
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Assembly in FASTA format
 *
 * @input database_dir
 * Optional path to database directory
 *
 * @input model_fp
 * Optional path to prodigal model file
 *
 * @output gbk      GenBank file of cap locus (optional)
 * @output svg      SVG visualization of cap locus (optional)
 * @output tsv      Tab-delimited hicap results
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process HICAP {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>
    database_dir   : Path?
    model_fp       : Path?

    output:
    gbk      = tuple(meta, files("*.gbk", optional: true))
    svg      = tuple(meta, files("*.svg", optional: true))
    tsv      = tuple(meta, file("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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
    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fasta.toList()[0].getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.toList()[0].getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    hicap \\
        --query_fp ${fasta_name} \\
        ${database_args} \\
        ${model_args} \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        --debug \\
        -o ./

    if [ ! -f ${prefix}.tsv ]; then
        echo "isolate<TAB>predicted_serotype<TAB>attributes<TAB>genes_identified<TAB>locus_location<TAB>region_I_genes<TAB>region_II_genes<TAB>region_III_genes<TAB>IS1016_hits" | sed 's/<TAB>/\t/g' > ${prefix}.tsv
        echo "${prefix}<TAB>cap_not_found<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-" | sed 's/<TAB>/\t/g' >> ${prefix}.tsv
    else
        sed -i 's/#isolate/isolate/' ${prefix}.tsv
    fi

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
