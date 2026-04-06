/**
 * Predict *Escherichia coli* serotype (O and H antigens).
 *
 * Uses [ECTyper](https://github.com/phac-nml/ectyper) to identify the O-antigen (lipopolysaccharide)
 * and H-antigen (flagella) genes in *E. coli* genome assemblies via BLAST. It provides a
 * standardized serotype prediction (e.g., O157:H7).
 *
 * @status stable
 * @keywords bacteria, escherichia coli, serotype, o-antigen, h-antigen, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation ectyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, txt, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited E. coli serotype predictions (O and H antigens)
 * - `txt`: BLAST allele details for serotype determination
 */
nextflow.preview.types = true

process ECTYPER {
    tag "${prefix}"
    label 'process_medium'

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
        txt: file("${prefix}.blast_alleles.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.blast_alleles.txt")
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

    ectyper \\
        ${task.ext.args} \\
        --cores ${task.cpus} \\
        --output ./ \\
        --input ${fna_name}

    # Cleanup
    mv output.tsv ${prefix}.tsv
    mv blastn_output_alleles.txt ${prefix}.blast_alleles.txt
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ectyper: \$(echo \$(ectyper --version 2>&1)  | sed 's/.*ectyper //; s/ .*\$//')
    END_VERSIONS
    """
}
