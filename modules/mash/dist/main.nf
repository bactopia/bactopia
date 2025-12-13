/**
 * Calculate genomic distances using MinHash sketches.
 *
 * Uses [Mash](https://github.com/marbl/Mash) to compute the distance between query sequences
 * and a reference database. It uses MinHash sketches to rapidly estimate the Jaccard index,
 * providing a fast approximation of Average Nucleotide Identity (ANI).
 *
 * @status stable
 * @keywords mash, distance, minhash, ani, comparison, taxonomy
 * @tags complexity:moderate input-type:multiple output-type:single features:conditional-logic
 * @citation mash
 *
 * @input tuple(meta, query)
 * - `meta`: Groovy Map containing sample information
 * - `query`: FASTA, FASTQ, or Mash sketch file to be queried
 *
 * @input reference
 * The reference file (FASTA, FASTQ, or Mash sketch) to compare against
 *
 * @output dist      A tab-delimited summary of the Mash distances and p-values
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MASH_DIST {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, query) : Tuple<Map, Set<Path>>
    reference      : Path

    output:
    dist     = tuple(meta, files("*.txt"))
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
    def is_compressed = reference.getName().endsWith(".xz") ? true : false
    def reference_name = reference.getName().replace(".xz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        xz -c -d ${reference} > ${reference_name}
    fi

    echo "reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes" | sed 's/<TAB>/\t/g' > ${prefix}-dist.txt
    mash \\
        dist \\
        -p ${task.cpus} \\
        ${task.ext.args} \\
        ${reference_name} \\
        ${query} | sed 's/.fna.gz//g' | sort -rn -k5,5 -t\$'\t' >> ${prefix}-dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
    END_VERSIONS
    """
}
