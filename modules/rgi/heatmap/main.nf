/**
 * Creates heatmaps from RGI JSON output files.
 *
 * This process executes rgi_heatmap to perform analysis
 *
 * @status stable
 * @keywords resistance, antimicrobial resistance, CARD, RGI, heatmap, visualization
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation rgi_heatmap
 *
 * @input tuple(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: RGI JSON output files
 *
 * @output heatmap  Heatmap files in various formats (CSV, EPS, PNG)
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
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
