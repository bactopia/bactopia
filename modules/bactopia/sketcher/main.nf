/**
 * ${$MODULE_DESCRIPTION}
.
 *
 * This process executes sketcher to perform analysis
 *
 * @status stable
 * @keywords ${MODULE_KEYWORDS}
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent
 * @citation sketcher
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Input file
 *
 * @input mash_db
 * Path parameter for mash_db
 *
 * @input sourmash_db
 * Path parameter for sourmash_db
 *
 * @output sig      Sig
 * @output msh      Msh
 * @output mash     Mash
 * @output sourmash Sourmash
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process SKETCHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>
    mash_db        : Path
    sourmash_db    : Path

    output:
    sig      = tuple(meta, file("${prefix}.sig"))
    msh      = tuple(meta, files("${prefix}-k*.msh"))
    mash     = tuple(meta, file("${prefix}-mash-refseq88-k21.txt"))
    sourmash = tuple(meta, file("${prefix}-sourmash-gtdb-rs207-k31.txt"))
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
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = mash_db.getName().endsWith(".xz") ? true : false
    def mash_name = mash_db.getName().replace(".xz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        xz -c -d ${mash_db} > ${mash_name}
    fi

    gzip -cd ${fasta} | mash sketch -o ${prefix}-k21 -k 21 ${task.ext.args} -I ${prefix} -
    gzip -cd ${fasta} | mash sketch -o ${prefix}-k31 -k 31 ${task.ext.args} -I ${prefix} -
    sourmash sketch dna ${task.ext.args2} --merge ${prefix} -o ${prefix}.sig ${fasta}

    # Mash Screen
    echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/\t/g' > ${prefix}-mash-refseq88-k21.txt
    gzip -cd ${fasta} | mash screen ${task.ext.args3} -p ${task.cpus} ${mash_name} - | sort -gr >> ${prefix}-mash-refseq88-k21.txt

    # Sourmash classify
    sourmash lca classify --query ${prefix}.sig --db ${sourmash_db} > ${prefix}-sourmash-gtdb-rs207-k31.txt

    # Cleanup
    rm -rf ${mash_name}

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/.*sourmash //;')
    END_VERSIONS
    """
}
