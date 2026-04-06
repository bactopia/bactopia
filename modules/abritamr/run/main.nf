/**
 * Detect antimicrobial resistance and virulence genes.
 *
 * Uses [abriTAMR](https://github.com/MDU-PHL/abritamr), a NATA (National Association of
 * Testing Authorities) accredited pipeline, to report the presence of reportable AMR
 * genes. It acts as a wrapper for AMRFinderPlus, formatted for clinical reporting standards
 * used in Victoria, Australia.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, antimicrobial resistance, nata, amrfinderplus
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output,compression,conditional-logic
 * @citation abritamr
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, summary?, matches, partials, virulence, amrfinder, results, logs, nf_logs, versions)
 * - `summary?`: Tab-delimited NATA-accredited AMR report summary
 * - `matches`: Tab-delimited list of matched AMR genes
 * - `partials`: Tab-delimited list of partially matched AMR genes
 * - `virulence`: Tab-delimited list of detected virulence genes
 * - `amrfinder`: Raw AMRFinderPlus output
 */
nextflow.preview.types = true

process ABRITAMR_RUN {
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
        summary: file("${prefix}.abritamr.tsv", optional: true),
        matches: file("${prefix}.summary_matches.tsv"),
        partials: file("${prefix}.summary_partials.tsv"),
        virulence: file("${prefix}.summary_virulence.tsv"),
        amrfinder: file("${prefix}.amrfinder.out"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.abritamr.tsv", optional: true),
            files("${prefix}.summary_matches.tsv"),
            files("${prefix}.summary_partials.tsv"),
            files("${prefix}.summary_virulence.tsv"),
            files("${prefix}.amrfinder.out")
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
        gzip -c -d ${fna} > ./${fna_name}
    fi

    abritamr run \\
        --contigs ${fna_name} \\
        --prefix ${prefix} \\
        ${task.ext.args} \\
        --jobs ${task.cpus}

    # Rename output files to prevent name collisions
    mv ${prefix}/summary_matches.txt ./${prefix}.summary_matches.tsv
    mv ${prefix}/summary_partials.txt ./${prefix}.summary_partials.tsv
    mv ${prefix}/summary_virulence.txt ./${prefix}.summary_virulence.tsv
    mv ${prefix}/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f ${prefix}/abritamr.txt ]; then
        # This file is not always present
        mv ${prefix}/abritamr.txt ./${prefix}.abritamr.tsv
    fi

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi
    rm -rf ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
