/**
 * Use RGI (Resistance Gene Identifier) to predict resistome(s) from protein or nucleotide data.
 *
 * This process executes rgi_main to perform analysis
 *
 * @status stable
 * @keywords resistance, antimicrobial resistance, CARD, RGI
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation rgi_main
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: FASTA file containing nucleotide or protein sequences
 *
 * @output json     RGI results in JSON format
 * @output tsv      RGI results in tab-separated format
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process RGI_MAIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    json     = tuple(meta, files("*.json", optional: true))
    tsv      = tuple(meta, files("*.txt"))
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
    rgi \\
        main \\
        ${task.ext.args} \\
        --clean \\
        --data wgs \\
        --num_threads ${task.cpus} \\
        --output_file ${prefix} \\
        --input_sequence ${fasta}

    # Remove empty json files
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
