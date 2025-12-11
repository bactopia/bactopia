/**
 * phyloFlash is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an illumina (meta)genomic dataset..
 *
 * This process executes phyloflash to perform analysis
 *
 * @status stable
 * @keywords metagenomics, illumina datasets, phylogenetic composition
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation phyloflash
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Channel containing single or paired-end reads
 *
 * @input _silva_db
 * Path parameter for _silva_db
 *
 * @input _univec_db
 * Path parameter for _univec_db
 *
 * @output supplemental Supplemental
 * @output aln          Alignment file
 * @output summary      JSON summary file
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
nextflow.preview.types = true

process PHYLOFLASH {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>
    _silva_db       : Path
    _univec_db      : Path

    output:
    supplemental = tuple(meta, files("${prefix}/*"))
    aln          = file("${prefix}/${prefix}.toalign.fasta", optional: true)
    summary      = file("${prefix}/${prefix}.phyloFlash.json", optional: true)
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

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
    def read_opts = meta.single_end ? "-read1 ${reads[0]}" : "-read1 ${reads[0]} -read2 ${reads[1]}"
    """
    mkdir ${prefix}
    phyloFlash.pl \\
        ${task.ext.args} \\
        ${read_opts} \\
        -lib ${prefix} \\
        -dbhome . \\
        -CPUs ${task.cpus}

    jsonify-phyloflash.py ${prefix}.phyloFlash > ${prefix}.phyloFlash.json
    mv ${prefix}.* ${prefix}


    if phyloflash-summary.py ${prefix}/ | grep -q -c "WARNING: Multiple SSUs were assembled by SPAdes"; then
        MULTI="1"
    fi

    if [ "${task.ext.allow_multiple_16s}" == "true" ]; then
        MULTI="0"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloFlash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
