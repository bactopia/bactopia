/**
 * Identify genes by local assembly of reads.
 *
 * Uses [ARIBA](https://github.com/sanger-pathogens/ariba) (Antimicrobial Resistance Identification
 * By Assembly) to detect AMR and virulence genes by creating local assemblies from paired-end reads.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords fastq, local assembly, antimicrobial resistance, virulence, ariba
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output,compression
 * @citation ariba
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by ARIBA)
 * - `lr`: Long reads (not supported by ARIBA)
 *
 * @input db
 * An [ARIBA](https://github.com/sanger-pathogens/ariba) prepared database
 *
 * @output report       A tab-delimited report of the [ARIBA](https://github.com/sanger-pathogens/ariba) analysis
 * @output summary      A comma-delimited summary of the report created using `ariba summary`
 * @output supplemental Supplemental output files from [ARIBA](https://github.com/sanger-pathogens/ariba)
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process ARIBA_RUN {
    tag "${prefix} - ${db_name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>
    db                      : Path

    output:
    report       = tuple(meta, file("${prefix}-report.tsv"))
    summary      = tuple(meta, file("${prefix}-summary.csv"))
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, files("versions.yml"))

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

    // Determine read type from explicit slots (ARIBA requires paired-end reads)
    has_r1 = r1 != null
    has_r2 = r2 != null

    def db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    mv ${db_name} ${db_name}db
    ariba \\
        run \\
        ${db_name}db/ \\
        ${r1} ${r2} \\
        ${db_name} \\
        ${task.ext.args} \\
        --threads ${task.cpus}

    ariba \\
        summary \\
        ${db_name}/summary \\
        ${db_name}/report.tsv \\
        --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \\
        --col_filter n \\
        --row_filter n

    # Rename to avoid naming collisions
    mv ${db_name}/report.tsv ./${prefix}-report.tsv
    mv ${db_name}/summary.csv ./${prefix}-summary.csv
    mv ${db_name}/ supplemental/

    # Cleanup
    rm -rf ${db_name}db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
