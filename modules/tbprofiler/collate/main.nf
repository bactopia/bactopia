/**
 * Collate TB-Profiler results from multiple samples.
 *
 * Uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to aggregate profiling results
 * from multiple samples into summary tables and files suitable for phylogenetic visualization.
 *
 * @status stable
 * @keywords tuberculosis, mycobacterium, drug resistance, collate, summary
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation tbprofiler
 *
 * @input record(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: List of TB-Profiler JSON output files
 *
 * @output record(meta, csv, variants_csv, variants_txt, itol, results, logs, nf_logs, versions)
 * - `csv`: Main collated results in CSV format
 * - `variants_csv`: Collated variants in CSV format
 * - `variants_txt`: Collated variants in text format
 * - `itol`: iTOL formatted files for visualization
 */
nextflow.preview.types = true

process TBPROFILER_COLLATE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, json: Set<Path>): Record

    stage:
    stageAs 'results-tmp/*', json

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        csv: file("tbprofiler.csv"),
        variants_csv: file("tbprofiler.variants.csv"),
        variants_txt: file("tbprofiler.variants.txt"),
        itol: files("*.itol.*.txt", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("tbprofiler.csv"),
            files("tbprofiler.variants.csv"),
            files("tbprofiler.variants.txt"),
            files("*.itol.*.txt", optional: true)
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

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler collate --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
