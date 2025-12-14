/**
 * Capsular typing of Streptococcus pneumoniae from Illumina reads.
 *
 * Uses [PneumoCaT](https://github.com/ukhsa-collaboration/PneumoCaT) (Pneumococcal Capsular Typing)
 * to assign capsular types to *Streptococcus pneumoniae* using a two-step approach: first matching
 * reads to a global database, then using a mapped-based approach for specific serogroup differentiation.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords pneumocat, streptococcus pneumoniae, capsular typing, serotyping
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation pneumocat
 *
 * @note
 * Negative results will cause non-0 exit codes from PneumoCaT
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by PneumoCaT)
 * - `lr`: Long reads (not supported by PneumoCaT)
 *
 * @output xml      The PneumoCaT result files in XML format
 * @output txt      A file containing the coverage information across the genes
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process PNEUMOCAT {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>

    output:
    xml      = tuple(meta, files("*.xml", optional: true))
    txt      = tuple(meta, files("*.coverage_summary.txt", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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

    // Determine read type from explicit slots (PneumoCaT requires paired-end reads)
    has_r1 = r1 != null
    has_r2 = r2 != null
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
