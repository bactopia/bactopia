/**
 * Pan-genome wide association studies.
 *
 * Uses [Scoary](https://github.com/AdmiralenOla/Scoary) to score the components of the pan-genome
 * for associations to specified traits (phenotypes). It is designed to work with the gene
 * presence/absence output from Roary.
 *
 * @status stable
 * @keywords scoary, pangenome, gwas, association, bacteria, roary
 * @tags complexity:simple input-type:multiple output-type:single features:conditional-logic
 * @citation scoary
 *
 * @input tuple(meta, genes)
 * - `meta`: Groovy Map containing sample information
 * - `genes`: Gene presence/absence CSV file (typically from Roary)
 *
 * @input traits
 * CSV file containing trait information for the samples
 *
 * @output csv      Scoary results files (*.csv)
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
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
    results  = tuple(meta, files("*.results.csv"))
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
