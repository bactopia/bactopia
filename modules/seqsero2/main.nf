/**
 * Salmonella serotype prediction from genome sequencing data.
 *
 * Uses [SeqSero2](https://github.com/denglab/SeqSero2) to predict *Salmonella* serotypes from
 * raw sequencing reads or genome assemblies using specific O-antigen and H-antigen markers.
 *
 * @status stable
 * @keywords salmonella, serotype, prediction, seqsero2, antigen
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation seqsero2
 *
 * @input tuple(meta, seqs)
 * - `meta`: Groovy Map containing sample information
 * - `seqs`: FASTQ reads or Assembled contigs
 *
 * @output tsv      SeqSero2 results in TSV format
 * @output txt      SeqSero2 results in text format
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SEQSERO2 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, seqs) : Tuple<Map, Set<Path>>

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    txt      = tuple(meta, files("${prefix}.txt"))
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
    def is_compressed_fna = seqs.toList()[0].getName().endsWith("fna.gz") ? true : false
    def seq_name = is_compressed_fna ? seqs.toList()[0].getName().replace(".gz", "") : "${seqs}"
    """
    if [ "${is_compressed_fna}" == "true" ]; then
        gzip -c -d ${seqs} > ${seq_name}
    fi

    SeqSero2_package.py \\
        ${task.ext.args} \\
        -d supplemental/ \\
        -n ${prefix} \\
        -p ${task.cpus} \\
        -t 4 \\
        -i ${seq_name}

    mv supplemental/SeqSero_log.txt ./${prefix}.log
    mv supplemental/SeqSero_result.tsv ./${prefix}.tsv
    mv supplemental/SeqSero_result.txt ./${prefix}.txt

    # Cleanup
    rm -rf supplemental/ ${seq_name} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}
