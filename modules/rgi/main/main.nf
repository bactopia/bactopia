/**
 * Predict antibiotic resistance from assemblies.
 *
 * Uses [RGI](https://github.com/arpcard/rgi) (Resistance Gene Identifier) to predict
 * resistomes from protein or nucleotide data based on homology and SNP models using
 * the Comprehensive Antibiotic Resistance Database (CARD).
 *
 * @status stable
 * @keywords resistance, antimicrobial resistance, card, rgi, amr
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation rgi_main
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output json     RGI results in JSON format
 * @output tsv      RGI results in tab-separated format
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process RGI_MAIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    json     = tuple(meta, file("${prefix}.json", optional: true))
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
    """
    rgi \\
        main \\
        ${task.ext.args} \\
        --clean \\
        --data wgs \\
        --num_threads ${task.cpus} \\
        --output_file ${prefix} \\
        --input_sequence ${assembly}

    # Remove empty json files
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi
    mv ${prefix}.txt ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
