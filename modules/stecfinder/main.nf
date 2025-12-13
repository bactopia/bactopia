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
 * @input tuple(meta, assembly, reads)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 * - `reads`: FASTQ reads (Illumina or Nanopore)
 *
 * @output tsv      TSV file with STEC gene markers results
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process STECFINDER {
    tag "${prefix}"
    label 'process_low'
    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly, reads) : Tuple<Map, Set<Path>, Set<Path>>

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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
    def is_compressed = (assembly.toList()[0].getName().endsWith(".gz") ? true : false) && !task.ext.stecfinder_use_reads ? true : false
    def seq_name = is_compressed ? assembly.toList()[0].getName().replace(".gz", "") : reads.join(" ")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${seq_name}
    fi

    stecfinder \\
        -i ${seq_name} \\
        ${task.ext.args} \\
        -t ${task.cpus} > ${prefix}.tsv

    # Cleanup
    rm -rf ${seq_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stecfinder: \$(echo \$(stecfinder --version 2>&1) | sed 's/^.*STECFinder version: //;' )
    END_VERSIONS
    """
}
