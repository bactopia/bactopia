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
 * @citation rgi_heatmap
 *
 * @input tuple(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: List of RGI results in JSON format
 *
 * @output heatmap  Heatmap files in various formats (CSV, EPS, PNG)
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process RGI_HEATMAP {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, json) : Tuple<Map, Set<Path>>

    stage:
    stageAs 'json/*', json

    output:
    heatmap  = tuple(meta, files("*.{csv,eps,png}", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "merged-results/"
    meta.logs_dir = "merged-results/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    NUM_SAMPLES=\$(ls json/ | wc -l)
    if [[ "\${NUM_SAMPLES}" -gt 1 ]]; then
        rgi \\
            heatmap \\
            ${task.ext.args} \\
            --output ${prefix} \\
            --input json/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
