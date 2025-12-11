/**
 * A tool to estimate bacterial species abundance.
 *
 * This process executes midas_species to perform analysis
 *
 * @status stable
 * @keywords bacteria, metagenomic, abundance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent
 * @citation midas_species
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Reads in FASTQ format
 *
 * @input db
 * A database formatted for MIDAS
 *
 * @output tsv        Tab-delimited summary of MIDAS results
 * @output abundances Species abundance profile from MIDAS
 * @output logs       Optional tool execution logs
 * @output nf_logs    Nextflow execution logs
 * @output versions   Software version information (YAML format)
 */
nextflow.preview.types = true

process MIDAS_SPECIES {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>
    db             : Path

    output:
    tsv        = tuple(meta, file("${prefix}.midas.tsv"))
    abundances = tuple(meta, files("*.abundances.txt"))
    logs       = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs   = tuple(meta, files(".command.*"))
    versions   = tuple(meta, file("versions.yml"))

    script:
    def VERSION = '1.3.2'
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
    def read_opts = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        MIDAS_DB=\$(find database/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    else
        MIDAS_DB=\$(find ${db}/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    fi

    run_midas.py \\
        species \\
        results \\
        ${read_opts} \\
        ${task.ext.args} \\
        -d \${MIDAS_DB} \\
        -t ${task.cpus}

    mv results/species/species_profile.txt ${prefix}.midas.abundances.txt
    midas-summary.py ${prefix} ${prefix}.midas.abundances.txt

    # Cleanup
    rm -rf results/
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: ${VERSION}
    END_VERSIONS
    """
}
