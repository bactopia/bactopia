/**
 * Predict O-antigen serogroup of Pseudomonas aeruginosa isolates.
 *
 * Uses [Pasty](https://github.com/rpetit3/pasty) (in silico serogrouping of *Pseudomonas aeruginosa* isolates)
 * to predict the O-antigen serogroup by searching the genome assembly for specific serogroup-associated
 * genes within the O-antigen locus.
 *
 * @status stable
 * @keywords bacteria, pseudomonas aeruginosa, serogroup, o-antigen, typing, blast
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation pasty
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv      A tab-delimited summary file with the predicted O-antigen serogroup
 * @output blast    A tab-delimited file of all raw BLAST hits used for the prediction
 * @output details  A tab-delimited file with detailed gene hits for each serogroup tested
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process PASTY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    blast    = tuple(meta, file("${prefix}.blastn.tsv"))
    details  = tuple(meta, file("${prefix}.details.tsv"))
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
    pasty \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${assembly}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/.*pasty, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pasty --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
