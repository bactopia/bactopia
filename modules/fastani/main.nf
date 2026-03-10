/**
 * Compute whole-genome Average Nucleotide Identity (ANI).
 *
 * Uses [FastANI](https://github.com/ParBLiSS/FastANI) to perform alignment-free computation
 * of ANI between the input query genomes and a reference genome. This is the standard method
 * for species definition (typically >95% ANI) and is much faster than traditional BLAST-based approaches.
 *
 * @status stable
 * @keywords fastani, ani, average nucleotide identity, taxonomy, genomic distance, comparison
 * @tags complexity:moderate input-type:multiple output-type:single features:conditional-logic
 * @citation fastani
 *
 * @input tuple(meta, query)
 * - `meta`: Groovy Map containing sample information
 * - `query`: One or more assembled contigs in FASTA format (Query genomes)
 *
 * @input reference
 * The reference genome assembly in FASTA format to compare against
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited summary of ANI scores, matched fragments, and total fragments
 */
nextflow.preview.types = true

process FASTANI {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, query: Set<Path>): Record
    reference      : Path

    stage:
    stageAs 'query-tmp/*', query

    output:
    record(
        meta: meta,
        tsv: file("${reference_name}.tsv"),
        results: [file("${reference_name}.tsv")],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    reference_fasta = reference.getName().replace(".gz", "")
    reference_name = reference_fasta.replace(".fna", "")
    prefix = task.ext.prefix ?: "${reference_name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "fastani/logs/${prefix}"
    meta.process_name = task.ext.process_name
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${reference_fasta}
    fi

    mkdir query
    cp -L query-tmp/* query/
    find query/ -name "*.gz" | xargs gunzip
    find query/ -name "*" -type f > query-list.txt

    fastANI \\
        --ql query-list.txt \\
        -r ${reference_fasta} \\
        -o fastani-result.tmp

    echo "query<TAB>reference<TAB>ani<TAB>mapped_fragments<TAB>total_fragments" | sed 's/<TAB>/\t/g' > ${reference_name}.tsv
    sed 's=^query/==' fastani-result.tmp >> ${reference_name}.tsv

    # Cleanup
    rm -rf ${reference_fasta} query/ query-list.txt fastani-result.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
