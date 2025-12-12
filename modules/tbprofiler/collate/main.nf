/**
 * Collate TB-Profiler results from multiple samples.
 *
 * Uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to aggregate profiling results
 * from multiple samples into summary tables and files suitable for phylogenetic visualization.
 *
 * @status stable
 * @keywords tuberculosis, mycobacterium, drug resistance, collate, summary
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation tbprofiler_collate
 *
 * @input tuple(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: List of TB-Profiler JSON output files
 *
 * @output csv          Main collated results in CSV format
 * @output variants_csv Collated variants in CSV format
 * @output variants_txt Collated variants in text format
 * @output itol         iTOL formatted files for visualization
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process TBPROFILER_COLLATE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, json) : Tuple<Map, Set<Path>>

    stage:
    stageAs 'results-tmp/*', json

    output:
    csv          = tuple(meta, file("tbprofiler.csv"))
    variants_csv = tuple(meta, file("tbprofiler.variants.csv"))
    variants_txt = tuple(meta, file("tbprofiler.variants.txt"))
    itol         = tuple(meta, files("*.itol.*.txt", optional: true))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = "merged-results"
    meta.logs_dir = "merged-results/logs/${meta.process_name}"
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    # Uncompress the JSON files
    mkdir results
    cp -L results-tmp/* results/
    find results/ -name "*.json.gz" | xargs gunzip

    tb-profiler \\
        collate \\
        ${task.ext.args} \\
        --db_dir database/ \\
        --format csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler collate --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
