/**
 * Pan-genome wide association studies.
 *
 * Uses [Scoary](https://github.com/AdmiralenOla/Scoary) to score the components of the pan-genome
 * for associations to specified traits (phenotypes). It is designed to work with the gene
 * presence/absence output from Roary.
 *
 * @status stable
 * @keywords scoary, pangenome, gwas, association, bacteria, roary
 * @tags complexity:simple input-type:multiple output-type:single features:conditional-logic
 * @citation scoary
 *
 * @input record(meta, csv)
 * - `meta`: Groovy Map containing sample information
 * - `csv`: Gene presence/absence CSV file (typically from Roary)
 *
 * @input traits
 * CSV file containing trait information for the samples
 *
 * @output record(meta, results, logs, nf_logs, versions)
 */
nextflow.preview.types = true

process SCOARY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, csv: Path): Record
    traits: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        // Generic fields (used for publishing)
        results: [
            files("*.results.csv")
        ],
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
    meta.process_name = task.ext.process_name
    meta.output_dir = "${meta.process_name}/"
    meta.logs_dir = "${meta.process_name}/logs"
    """
    scoary \\
        ${task.ext.args} \\
        --no-time \\
        --threads ${task.cpus} \\
        --traits ${traits} \\
        --genes ${csv}

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scoary: \$( scoary --version 2>&1 )
    END_VERSIONS
    """
}
