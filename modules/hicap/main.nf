/**
 * Predict *Haemophilus influenzae* capsule serotype.
 *
 * Uses [hicap](https://github.com/scwatts/hicap) to identify the capsule locus in *H. influenzae*
 * genome assemblies. It predicts the serotype (a, b, c, d, e, f, or Non-Typeable/NTHi) and
 * can optionally generate visualizations of the locus structure.
 *
 * @status stable
 * @keywords bacteria, haemophilus influenzae, serotype, capsule, typing, nthi
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation hicap
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input database_dir
 * Optional path to a custom hicap reference database directory
 *
 * @input model_fp
 * Optional path to a custom Prodigal training model file
 *
 * @output gbk       GenBank file containing the annotated capsule locus region (optional)
 * @output svg       SVG visualization of the capsule locus gene arrangement (optional)
 * @output tsv       A tab-delimited summary of the predicted serotype and locus coverage
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process HICAP {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>
    database_dir   : Path?
    model_fp       : Path?

    output:
    gbk      = tuple(meta, files("*.gbk", optional: true))
    svg      = tuple(meta, files("*.svg", optional: true))
    tsv      = tuple(meta, files("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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
    def is_compressed = assembly.toList()[0].getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.toList()[0].getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    hicap \\
        --query_fp ${assembly_name} \\
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
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
