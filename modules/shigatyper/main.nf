/**
 * Shigella serotype from Illumina or Oxford Nanopore reads.
 *
 * Uses [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper) to determine the serotype
 * of *Shigella* isolates using Illumina paired-end reads or Oxford Nanopore long reads. It detects
 * serotype-specific genes and markers to provide a predicted serotype.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords shigella, serotype, typing, illumina, nanopore, reads
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation shigatyper
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output record(meta, tsv, hits, results, logs, nf_logs, versions)
 * - `tsv`: ShigaTyper serotype prediction results in TSV format
 * - `hits`: Detailed gene hits from ShigaTyper analysis
 */
nextflow.preview.types = true

process SHIGATYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        hits: files("${prefix}-hits.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}-hits.tsv", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null
    meta.single_end = has_se && !has_r1 && !has_r2

    if (has_lr) {
        """
        shigatyper \\
            ${task.ext.args} \\
            --SE ${lr} \\
            --ont \\
            --name ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
    else if (meta.single_end) {
        """
        shigatyper \\
            ${task.ext.args}  \\
            --SE ${se} \\
            --name ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
    else {
        """
        shigatyper \\
            ${task.ext.args}  \\
            --R1 ${r1} \\
            --R2 ${r2} \\
            --name ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}
