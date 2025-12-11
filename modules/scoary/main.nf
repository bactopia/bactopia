/**
 * Pan-genome wide association studies.
 *
 * This process executes scoary to perform analysis
 *
 * @status stable
 * @keywords pangenome, GWAS, association, Roary
 * @tags complexity:simple input-type:multiple output-type:single
 * @citation scoary
 *
 * @input tuple(meta, genes)
 * - `meta`: Groovy Map containing sample information
 * - `genes`: Gene presence/absence file from Roary
 *
 * @input traits
 * Path parameter for traits
 *
 * @output csv      Scoary results files
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process SCOARY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, genes) : Tuple<Map, Path>
    traits         : Path

    output:
    csv      = tuple(meta, files("*.csv"))
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
    meta.process_name = task.ext.process_name
    meta.output_dir = "${meta.process_name}/"
    meta.logs_dir = "${meta.process_name}/logs"
    """
    scoary \\
        ${task.ext.args} \\
        --no-time \\
        --threads ${task.cpus} \\
        --traits ${traits} \\
        --genes ${genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scoary: \$( scoary --version 2>&1 )
    END_VERSIONS
    """
}
