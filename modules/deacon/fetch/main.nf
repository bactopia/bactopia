/**
 * Fetch a pre-built deacon index for host read filtering.
 *
 * Uses [deacon](https://github.com/bede/deacon) to download a pre-built minimizer index
 * for filtering host reads from sequencing data. The default index is `panhuman-1`, which
 * covers human reference genomes for human read depletion.
 *
 * @status stable
 * @keywords host, decontamination, depletion, download, index, minimizer, deacon
 * @tags complexity:simple input-type:none output-type:single features:internet-access,resource-download,no-test
 * @citation deacon
 *
 * @note Internet Required
 * This process requires an active internet connection to fetch the pre-built index.
 *
 * @output record(db, logs)
 * - `db`: The pre-built deacon minimizer index file
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028
process DEACON_FETCH {
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("${prefix}/${task.ext.deacon_index_name}.idx"),
        logs: files("${prefix}/logs/*", optional: true)
    )

    script:
    prefix = task.ext.process_name
    """
    mkdir -p ${prefix}/logs

    deacon \\
        index \\
        fetch \\
        ${task.ext.deacon_index_name} \\
        > ${prefix}/${task.ext.deacon_index_name}.idx

    # Move outputs to tool specific folder
    cp .command.begin ${prefix}/logs/nf.command.begin
    cp .command.err ${prefix}/logs/nf.command.err
    cp .command.log ${prefix}/logs/nf.command.log
    cp .command.out ${prefix}/logs/nf.command.out
    cp .command.run ${prefix}/logs/nf.command.run
    cp .command.sh ${prefix}/logs/nf.command.sh
    cp .command.trace ${prefix}/logs/nf.command.trace

    cat <<-END_VERSIONS > ${prefix}/logs/versions.yml
    "${task.process}":
        deacon: \$( deacon --version | head -n1 | sed 's/deacon //' )
    END_VERSIONS
    """
}
