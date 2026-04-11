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
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input database_dir?
 * Path to a custom hicap reference database directory
 *
 * @input model_fp?
 * Path to a custom Prodigal training model file
 *
 * @output record(meta, gbff?, svg?, tsv, results, logs, nf_logs, versions)
 * - `gbff?`: GenBank file containing the annotated capsule locus region
 * - `svg?`: SVG visualization of the capsule locus gene arrangement
 * - `tsv`: Tab-delimited summary of the predicted serotype and locus coverage
 */
nextflow.preview.types = true

process HICAP {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )
    database_dir: Path?
    model_fp    : Path?

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        gbff: file("${prefix}.gbk", optional: true),
        svg: file("${prefix}.svg", optional: true),
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.gbk", optional: true),
            files("${prefix}.svg", optional: true)
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
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )

    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    hicap \\
        --query_fp ${fna_name} \\
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
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
