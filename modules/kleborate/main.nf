/**
 * Screening genomic assemblies of Klebsiella for clinically relevant determinants.
 *
 * This process executes kleborate to perform analysis
 *
 * @status stable
 * @keywords klebsiella, resistance, virulence, typing
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation kleborate
 *
 * @input tuple(meta, fastas)
 * - `meta`: Groovy Map containing sample information
 * - `fastas`: Assembly files in FASTA format
 *
 * @output txt      Tab-delimited Kleborate results
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process KLEBORATE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fastas) : Tuple<Map, Path>

    output:
    txt      = tuple(meta, files("*.txt"))
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
    """
    mkdir results/
    kleborate \\
        ${task.ext.args} \\
        --outdir results/ \\
        --assemblies ${fastas}

    # Rename output file to include the prefix name
    find results/ -name "*output.txt" -print0 | while read -d \$'\0' file; do mv "\$file" "${prefix}.txt"; done

    # Negative results will not have an output file
    if [ ! -f "${prefix}.txt" ]; then
        touch "${prefix}.txt"
    fi

    # cleanup
    rm -rf results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
