/**
 * Predict *Haemophilus parasuis* serotype.
 *
 * Uses [HPSuisSero](https://github.com/Abraham-L/HPSuisSero) to predict the serotype of
 * *Haemophilus parasuis* (syn. *Glaesserella parasuis*) assemblies. It detects specific
 * capsule loci markers to assign one of the known serovars.
 *
 * @status stable
 * @keywords bacteria, haemophilus parasuis, glaesserella parasuis, serotype, typing, capsule
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation hpsuissero
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv       A tab-delimited summary of the predicted serotype and gene hits
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process HPSUISSERO {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
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

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    HpsuisSero.sh \\
        -i ${assembly_name} \\
        -o ./ \\
        -s ${prefix} \\
        -x fasta \\
        -t ${task.cpus}

    # Cleanup
    mv ${prefix}_serotyping_res.tsv ./${prefix}.tsv
    rm -rf ${assembly_name} blast_res/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpsuissero: ${task.ext.version}
    END_VERSIONS
    """
}
