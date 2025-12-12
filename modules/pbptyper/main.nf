/**
 * Predict Penicillin Binding Protein (PBP) type of *Streptococcus pneumoniae* assemblies.
 *
 * Uses [PBPtyper](https://github.com/rpetit3/pbptyper) to detect variations in the three
 * key PBP genes (*pbp1a*, *pbp2b*, and *pbp2x*) in *S. pneumoniae*. Typing these genes is
 * essential for predicting reduced susceptibility or full resistance to penicillin and other
 * $\beta$-lactam antibiotics.
 *
 * @status stable
 * @keywords bacteria, streptococcus pneumoniae, penicillin, amr, resistance, pbp, typing
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation pbptyper
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv       A tab-delimited summary file with the predicted PBP type for each gene
 * @output blast     A tab-delimited file of the raw TBLASTN hits used for gene identification
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process PBPTYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    blast    = tuple(meta, files("*.tblastn.tsv"))
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
    pbptyper \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${assembly}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/.*pbptyper, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pbptyper --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
