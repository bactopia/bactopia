/**
 * Capsular typing of Streptococcus pneumoniae from Illumina reads.
 *
 * This process executes pneumocat to perform analysis
 *
 * @status stable
 * @keywords pneumocat, Streptococcus pneumoniae, capsular typing, serotyping
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation pneumocat
 *
 * @input tuple(meta, meta)
 * - `meta`: Groovy Map containing sample information
 * - `meta`: Groovy Map containing sample information
 *
 * @output xml      The pneumocat result files in xml format
 * @output txt      A file containing the coverage information acrosss the genes
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process PNEUMOCAT {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, _reads) : Tuple<Map, Path>

    output:
    xml      = tuple(meta, files("*.xml", optional: true))
    txt      = tuple(meta, files("*.coverage_summary.txt", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    def VERSION = '1.2.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
    PneumoCaT.py \\
        --input_directory ./ \\
        --threads ${task.cpus} \\
        --output_dir ./

    # clean up
    rm -rf *.bam *.bai ComponentComplete.txt

    # PneumoCAT uses first match in a glob, so moves between R1 and R2
    if [ -f ${prefix}_R1.results.xml ]; then
        mv ${prefix}_R1.results.xml ${prefix}.results.xml
    else
        mv ${prefix}_R2.results.xml ${prefix}.results.xml
    fi
    mv logs/* ./
    mv coverage_summary.txt ${prefix}.coverage_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: ${VERSION}
    END_VERSIONS
    """
}
