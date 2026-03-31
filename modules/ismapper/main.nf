/**
 * Identify insertion sites and orientation of mobile genetic elements.
 *
 * Uses [ISMapper](https://github.com/jhawkey/IS_mapper) to identify the position and orientation
 * of specific insertion sequences (IS) in a bacterial genome. It works by mapping paired-end reads
 * to a library of IS queries and a reference genome to determine where the IS elements are located
 * relative to the reference coordinates.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2) where each read slot is Path
 *
 * @status stable
 * @keywords bacteria, mobile elements, insertion sequences, mapping, structural variation, ismapper
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation ismapper
 *
 * @input record(meta, r1, r2)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 *
 * @input reference
 * Reference genome in GenBank format (*.gbk) to map insertion sites against
 *
 * @input query
 * FASTA file containing the insertion sequences to search for
 *
 * @output record(meta, results, logs, nf_logs, versions)
 *
 * @results supplemental
 * - `*.bed.gz`: Insertion site coordinates in BED format
 * - `*.fastq.gz`: Supporting reads mapped to insertion sites
 */
nextflow.preview.types = true

process ISMAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, r1: Path, r2: Path): Record
    reference: Path
    query    : Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        // Generic fields (used for publishing)
        results: [
            files("supplemental/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    def query_name = query.getName().replace(".gz", "")
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${query_name}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${query_name}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    def ref_compressed = reference.getName().endsWith(".gz") ? true : false
    def reference_name = reference.getName().replace(".gz", "")
    def query_compressed = query.getName().endsWith(".gz") ? true : false
    """
    if [ "${ref_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${reference_name}
    fi
    if [ "${query_compressed}" == "true" ]; then
        gzip -c -d ${query} > ${query_name}
    fi

    ismap \\
        ${task.ext.args} \\
        --t ${task.cpus} \\
        --output_dir ${prefix} \\
        --queries ${query_name} \\
        --log ${prefix} \\
        --reference ${reference_name} \\
        --reads ${r1} ${r2}

    # Reorganize output files
    mkdir supplemental
    mv ${prefix}/*/* supplemental/

    # Cleanup and compress FASTQ and BED files
    rm -rf ${reference_name} ${query_name} ${prefix}/
    find supplemental/ -name "*.fastq" | xargs -I {} gzip {}
    find supplemental/ -name "*.bed" | xargs -I {} gzip {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ismapper: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
