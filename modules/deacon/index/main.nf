/**
 * Build a deacon minimizer index from a FASTA reference genome.
 *
 * Uses [deacon](https://github.com/bede/deacon) to build a minimizer index from a reference
 * genome in FASTA format. The resulting index is used by the deacon filter module to remove
 * host contamination from sequencing reads via SIMD-accelerated minimizer comparison.
 *
 * @status stable
 * @keywords host, decontamination, depletion, index, minimizer, reference, deacon
 * @tags complexity:simple input-type:single output-type:single features:no-test
 * @citation deacon
 *
 * @note Build Once
 * This process builds a minimizer index from a reference FASTA. The index is cached
 * via storeDir and only needs to be built once per reference genome.
 *
 * @input reference
 * Reference genome in FASTA format to build the minimizer index from
 *
 * @output record(db, logs)
 * - `db`: The built deacon minimizer index file
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028
process DEACON_INDEX {
    tag "deacon-index"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    reference: Path

    output:
    record(
        db: file("${prefix}/deacon-index.idx"),
        logs: files("${prefix}/logs/*", optional: true)
    )

    script:
    prefix = task.ext.process_name
    """
    mkdir -p ${prefix}/logs

    deacon \\
        index \\
        build \\
        --threads ${task.cpus} \\
        ${task.ext.args} \\
        ${reference} \\
        > ${prefix}/deacon-index.idx

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
