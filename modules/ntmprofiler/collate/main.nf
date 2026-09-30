/**
 * Collate NTM-Profiler results from multiple samples.
 *
 * Uses [NTM-Profiler](https://github.com/jodyphelan/NTM-Profiler) to merge
 * species, lineage, resistance, and variant results into run-level tables.
 *
 * @status stable
 * @keywords mycobacterium, ntm, species identification, antimicrobial resistance, summary
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent,aggregation
 * @citation ntmprofiler
 *
 * @input record(meta, json)
 * - `meta`: Groovy Record containing sample information
 * - `json`: NTM-Profiler JSON result files
 *
 * @input db
 * Directory or compressed tarball containing an NTM-Profiler database
 *
 * @output record(meta, csv, variants_csv, snp_dist_json?, results, logs, nf_logs, versions)
 * - `csv`: Collated species, lineage, and resistance results
 * - `variants_csv`: Collated variant calls
 * - `snp_dist_json?`: SNP-distance graph when a distance database is present
 */
nextflow.enable.types = true

process NTMPROFILER_COLLATE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record(meta: Record, json: Set<Path>)
    db: Path

    stage:
    stageAs json, 'staging/json/*'

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        csv: file("ntmprofiler.csv"),
        variants_csv: file("ntmprofiler.csv.variants.csv"),
        snp_dist_json: file("ntmprofiler.csv.snp-dist.json", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("ntmprofiler.csv"),
            files("ntmprofiler.csv.variants.csv"),
            files("ntmprofiler.csv.snp-dist.json", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"
    is_tarball = db.getName().endsWith(".tar.gz")
    db_setup = is_tarball ? "mkdir -p database && tar -xzf ${db} -C database" : ""
    db_arg = is_tarball ? "--db_dir database" : "--db_dir ${db}"
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        process_name: task.ext.process_name,
        output_dir: "merged-results",
        logs_dir: "merged-results/logs/${task.ext.process_name}"
    )
    """
    mkdir -p results
    ${db_setup}
    cp -L staging/json/* results/
    find results -name '*.json.gz' -exec gunzip {} \\;

    ntm-profiler \\
        collate \\
        ${task.ext.args} \\
        ${db_arg} \\
        --dir results \\
        --outfile ntmprofiler.csv \\
        --format csv

    # Cleanup
    rm -rf database results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ntm-profiler: \$(ntm-profiler --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}
