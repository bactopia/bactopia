/**
 * Determine the phylogroup of Escherichia coli isolates.
 *
 * Uses [ClermonTyping](https://github.com/A-BN/ClermonTyping) to perform in silico PCR
 * detection of specific marker genes (arpA, chuA, yjaA, TspE4.C2). This assigns the
 * isolate to one of the main *E. coli* phylogroups (A, B1, B2, C, D, E, F, G, or Cryptic).
 *
 * @status stable
 * @keywords bacteria, escherichia coli, typing, phylotyping, pcr, phylogroup
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation clermontyping
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited E. coli phylogroup assignment with detected marker genes
 *
 * @results supplemental
 * - `*.blast.xml`: Raw BLAST XML output for marker gene detection
 * - `*.mash.tsv`: Mash screen results for phylogroup confirmation
 * - `*.html`: Interactive HTML report of the phylogroup assignment
 */
nextflow.preview.types = true

process CLERMONTYPING {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        fna: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("supplemental/*")
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
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    clermonTyping.sh \\
        --fasta ${fna_name} \\
        --name supplemental \\
        ${task.ext.args}

    # Remove temporary files and rename outputs
    rm supplemental/${fna_name} ${fna_name}
    rm -rf supplemental/db
    rm -rf supplemental/supplemental.R
    mv supplemental/${fna_name}.xml supplemental/${prefix}.blast.xml
    mv supplemental/${fna_name}_mash_screen.tab supplemental/${prefix}.mash.tsv
    mv supplemental/supplemental.html supplemental/${prefix}.html

    # add column names to phylogroups file
    echo "sample<TAB>detected_genes<TAB>quadruplex_genes<TAB>quadruplex_alleles<TAB>phylogroup<TAB>mash_group" | sed 's/<TAB>/\t/g' > ${prefix}.tsv
    cat supplemental/supplemental_phylogroups.txt >> ${prefix}.tsv

    # Cleanup
    rm supplemental/supplemental_phylogroups.txt
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clermontyping: \$(echo \$(clermonTyping.sh -v 2>&1) | sed 's/^.* version : //;s/ .*\$//')
    END_VERSIONS
    """
}
