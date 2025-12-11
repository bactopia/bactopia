/**
 * Shigella serotyping from Illumina or Oxford Nanopore reads.
 *
 * This process executes shigatyper to perform analysis
 *
 * @status stable
 * @keywords Shigella, serotype, typing, reads
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation shigatyper
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ file(s) (single-end, paired-end, or ONT)
 *
 * @output tsv      ShigaTyper results in TSV format
 * @output hits     Detailed hits from ShigaTyper
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process SHIGATYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    hits     = tuple(meta, file("${prefix}-hits.tsv", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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

    if (_meta.runtype == "ont") {
        """
        shigatyper \\
            ${task.ext.args} \\
            --SE ${reads} \\
            --ont \\
            --name ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
    else if (_meta.single_end) {
        """
        shigatyper \\
            ${task.ext.args}  \\
            --SE ${reads} \\
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
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}
