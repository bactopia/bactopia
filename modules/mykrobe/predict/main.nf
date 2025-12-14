/**
 * Predict Antimicrobial Resistance (AMR) for supported bacterial species.
 *
 * Uses [Mykrobe](https://github.com/mykrobe/mykrobe) to quickly predict resistance and susceptibility
 * based on short reads (FASTQ) or aligned reads (BAM). It maps k-mers from the input sequences
 * against a curated database of resistance markers for species like *M. tuberculosis* and *S. aureus*.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords amr, resistance, susceptibility, k-mer, fastq, bam, mykrobe
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent
 * @citation mykrobe
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input species
 * The target species for which to make the AMR prediction (e.g., "tb" or "staph")
 *
 * @output csv       AMR predictions in machine-readable CSV format
 * @output json      Detailed AMR prediction results in JSON format
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MYKROBE_PREDICT {
    tag "${meta.name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>
    species                 : String

    output:
    csv      = tuple(meta, file("${prefix}.csv"))
    json     = tuple(meta, file("${prefix}.json"))
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

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs and options for mykrobe
    def is_ont = has_lr ? "--ont" : ""
    def read_inputs = has_lr ? "${lr}" : (meta.single_end ? "${se}" : "${r1} ${r2}")
    """
    mykrobe \\
        predict \\
        ${task.ext.args} ${is_ont} \\
        --species ${species} \\
        --threads ${task.cpus} \\
        --sample ${prefix} \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq ${read_inputs}

    # Cleanup
    rm -rf mykrobe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
