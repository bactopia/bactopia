/**
 * Serotype of Shigatoxin producing E. coli using reads/assemblies.
 *
 * Uses [STECFinder](https://github.com/LanLab/STECFinder) to identify Shiga toxin-producing
 * *Escherichia coli* (STEC) serotypes and virulence factors from genome assemblies or sequencing reads.
 *
 * @status stable
 * @keywords stec, e. coli, virulence, serotype, typing
 * @tags complexity:moderate input-type:single output-type:single features:conditional-logic
 * @citation stecfinder
 *
 * @input record(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: TSV file with STEC gene markers results
 */
nextflow.preview.types = true

process STECFINDER {
    tag "${prefix}"
    label 'process_low'
    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record

    stage:
    stageAs 'staging/fna/*', fna
    stageAs 'staging/r1/*', r1
    stageAs 'staging/r2/*', r2
    stageAs 'staging/se/*', se
    stageAs 'staging/lr/*', lr

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv")
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // Determine input type and construct sequence input
    def reads = [r1, r2, se, lr].findAll{ r -> r != null }
    def use_assembly = !task.ext.stecfinder_use_reads
    def is_compressed = use_assembly && fna.getName().endsWith(".gz")
    def seq_name = is_compressed ? fna.getName().replace(".gz", "") : use_assembly ? fna.getName() : reads.join(" ")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${seq_name}
    fi

    stecfinder \\
        -i ${seq_name} \\
        ${task.ext.args} \\
        -t ${task.cpus} > ${prefix}.tsv

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${seq_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stecfinder: \$(echo \$(stecfinder --version 2>&1) | sed 's/^.*STECFinder version: //;' )
    END_VERSIONS
    """
}
