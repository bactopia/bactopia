/**
 * Download and prepare the NTM-Profiler database.
 *
 * Uses `ntm-profiler update_db` to fetch and index the default NTM database once
 * for reuse by all profiling and collation tasks in a workflow run.
 *
 * @status stable
 * @keywords mycobacterium, ntm, database, download, antimicrobial resistance
 * @tags complexity:moderate input-type:none output-type:single features:database-dependent,resource-download
 * @citation ntmprofiler
 *
 * @note Internet Required
 * This process requires internet access when the cached database is unavailable.
 *
 * @output record(db, logs, nf_logs, versions)
 * - `db`: Prepared NTM-Profiler database directory
 */
nextflow.enable.types = true

// bactopia-lint: ignore M017,M018,M023,M024,M025,M026
process NTMPROFILER_UPDATE {
    tag "update"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("ntmprofiler-db"),
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    """
    ntm-profiler update_db \\
        --db_dir ntmprofiler-db \\
        --dir update-work

    # Cleanup
    rm -rf update-work

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ntm-profiler: \$(ntm-profiler --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}
