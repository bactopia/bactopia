/**
 * Predict Shigella serotypes and differentiate Shigella/EIEC.
 *
 * Uses [ShigaPass](https://github.com/Munch-Lab/ShigaPass) to identify *Shigella* serotypes and
 * distinguish *Shigella* species from Enteroinvasive *Escherichia coli* (EIEC) using specific
 * genomic markers from assembled contigs.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, virulence, prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation shigapass
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, flex_tsv, results, logs, nf_logs, versions)
 * - `tsv`: ShigaPass summary results in TSV format
 * - `flex_tsv`: ShigaPass Flex summary results in TSV format
 */
nextflow.preview.types = true

process SHIGAPASS {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        flex_tsv: file("${prefix}_Flex_summary.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}_Flex_summary.tsv", optional: true),
            files("supplemental/*", optional: true)
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

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")

    // WARN: Version information needed for -p parameter, see below. Please update this string when bumping container versions.
    def SHIGAPASS_VERSION = "1.5.0"
    """
    # ShigaPass does not accept compressed FASTA files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    # Convert our genome path to a file with a path in it
    ls ${fna_name} > ${fna_name}_tmp.txt

    # Run ShigaPass
    ShigaPass.sh \\
        -l ${fna_name}_tmp.txt \\
        ${task.ext.args} \\
        -p "\$(dirname \$(which ShigaPass.sh))/../share/shigapass-${SHIGAPASS_VERSION}/db" \\
        -t ${task.cpus} \\
        -o ${prefix}

    # Remove the temporary file from above
    if [ "${is_compressed}" == "true" ]; then
        rm ${fna_name}
    fi
    rm ${fna_name}_tmp.txt

    # Convert to tab delimited and move to the pwd
    mkdir -p supplemental
    sed 's/;/\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv
    mv ${prefix}/ShigaPass_summary.csv supplemental/

    # Convert to tab delimited and move to the pwd
    if [ -f ${prefix}/ShigaPass_Flex_summary.csv ]; then
        sed 's/;/\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv
        mv ${prefix}/ShigaPass_Flex_summary.csv supplemental/
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """
}
