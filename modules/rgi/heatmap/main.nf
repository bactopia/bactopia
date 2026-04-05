/**
 * Create heatmaps of resistance gene presence/absence.
 *
 * Uses [RGI](https://github.com/arpcard/rgi) (Resistance Gene Identifier) to generate
 * heatmaps visualizing the presence or absence of antimicrobial resistance genes across
 * multiple samples based on RGI JSON results.
 *
 * @status stable
 * @keywords resistance, antimicrobial resistance, card, rgi, heatmap, visualization
 * @tags complexity:simple input-type:multiple output-type:multiple features:conditional-logic
 * @citation rgi
 *
 * @input record(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: List of RGI results in JSON format
 *
 * @output record(meta, heatmap, results, logs, nf_logs, versions)
 * - `heatmap`: Heatmap files in various formats (CSV, EPS, PNG)
 */
nextflow.preview.types = true

process RGI_HEATMAP {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, json: Set<Path>): Record

    stage:
    stageAs 'staging/json/*', json

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        heatmap: files("*.{csv,eps,png}", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("*.{csv,eps,png}", optional: true)
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
    meta.output_dir = "rgi-heatmap"
    meta.logs_dir = "rgi-heatmap/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    NUM_SAMPLES=\$(ls staging/json/ | wc -l)
    if [[ "\${NUM_SAMPLES}" -gt 1 ]]; then
        rgi \\
            heatmap \\
            ${task.ext.args} \\
            --output ${prefix} \\
            --input staging/json/
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
