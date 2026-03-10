/**
 * Reconstruct and type plasmids from a bacterial genome assembly.
 *
 * Uses [MobSuite's mob_recon](https://github.com/phac-nml/mob-suite) to reconstruct plasmids
 * by grouping relevant contigs. It then uses the Mob-typer component to classify the plasmids
 * based on replicon type, incompatibility group (Inc type), and predicted mobility.
 *
 * @status stable
 * @keywords bacteria, plasmid, reconstruction, mobtyper, replicon, contigs, assembly
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation mobsuite
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, chromosome, contig_report, txt, plasmids, results, logs, nf_logs, versions)
 * - `chromosome`: Chromosomal sequences separated from plasmid contigs (gzipped FASTA)
 * - `contig_report`: Tab-delimited report assigning each contig to chromosome or plasmid
 * - `txt`: MOB-typer results with replicon type, mobility, and incompatibility group (optional)
 * - `plasmids`: Reconstructed plasmid sequences in gzipped FASTA format (optional)
 */
nextflow.preview.types = true

process MOBSUITE_RECON {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        chromosome: file("${prefix}-chromosome.fasta.gz"),
        contig_report: file("${prefix}-contig_report.txt"),
        txt: file("${prefix}-mobtyper.txt", optional: true),
        plasmids: files("plasmid_*.fasta.gz", optional: true),
        // Generic fields (used for publishing)
        results: files("${prefix}-chromosome.fasta.gz") + files("${prefix}-contig_report.txt") + files("${prefix}-mobtyper.txt", optional: true) + files("plasmid_*.fasta.gz", optional: true),
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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
    def fasta_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${fasta_name}
    fi

    mob_recon \\
        --infile ${fasta_name} \\
        ${task.ext.args} \\
        --num_threads ${task.cpus} \\
        --outdir supplemental \\
        --sample_id ${prefix}

    if [[ -f "supplemental/mobtyper_results.txt" ]]; then
        mv supplemental/mobtyper_results.txt ${prefix}-mobtyper.txt
    fi

    if [[ -f "supplemental/chromosome.fasta" ]]; then
        mv supplemental/chromosome.fasta ${prefix}-chromosome.fasta
        gzip ${prefix}-chromosome.fasta
    fi

    if [[ -f "supplemental/contig_report.txt" ]]; then
        mv supplemental/contig_report.txt ${prefix}-contig_report.txt
    fi

    # Cleanup
    gzip supplemental/*.fasta
    mv supplemental/*.fasta.gz ./
    rm -rf ${fasta_name} supplemental/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
