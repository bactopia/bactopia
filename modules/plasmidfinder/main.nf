/**
 * Identify plasmid replicon types in bacterial sequences and assemblies.
 *
 * Uses [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/)
 * to identify plasmid types (replicon typing) by querying the genome assembly against a database
 * of plasmid sequences. This is a crucial step for understanding the mobility of resistance and virulence genes.
 *
 * @status stable
 * @keywords plasmid, replicon, typing, identification, mobility, amr
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation plasmidfinder
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, json, txt, tsv, genome_seq, plasmid_seq, results, logs, nf_logs, versions)
 * - `json`: PlasmidFinder results in JSON format
 * - `txt`: PlasmidFinder results in text format
 * - `tsv`: Tab-delimited PlasmidFinder results with replicon typing information
 * - `genome_seq`: FASTA sequences of plasmid hits found in the genome (gzipped)
 * - `plasmid_seq`: Reference plasmid sequences matched (gzipped)
 */
nextflow.preview.types = true

process PLASMIDFINDER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        json: file("${prefix}.json"),
        txt: file("${prefix}.txt"),
        tsv: file("${prefix}.tsv"),
        genome_seq: file("${prefix}-hit_in_genome_seq.fsa.gz"),
        plasmid_seq: file("${prefix}-plasmid_seqs.fsa.gz"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.json"),
            files("${prefix}.txt"),
            files("${prefix}.tsv"),
            files("${prefix}-hit_in_genome_seq.fsa.gz"),
            files("${prefix}-plasmid_seqs.fsa.gz")
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

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '2.1.6'
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    plasmidfinder.py \\
        ${task.ext.args} \\
        -i ${fna_name} \\
        -o ./ \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv data.json ${prefix}.json
    mv results.txt ${prefix}.txt
    mv Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${prefix}-plasmid_seqs.fsa

    # Add sample name to TSV results
    head -n 1 results_tab.tsv | sed "s/^/Sample\t/" > ${prefix}.tsv
    tail -n +2 results_tab.tsv | sed "s/^/${prefix}\t/" >> ${prefix}.tsv

    # Cleanup
    gzip *.fsa
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi
    rm -rf results_tab.tsv tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: ${VERSION}
    END_VERSIONS
    """
}
