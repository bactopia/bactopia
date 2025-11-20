nextflow.preview.types = true

process SKETCHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    mash_db        : Path
    sourmash_db    : Path

    output:
    sig      = tuple(meta, file("${prefix}.sig"))
    msh      = tuple(meta, file("${prefix}-k*.msh"))
    mash     = tuple(meta, file("${prefix}-mash-refseq88-k21.txt"))
    sourmash = tuple(meta, file("${prefix}-sourmash-gtdb-rs207-k31.txt"))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_out   = tuple(meta, file(".command.out"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
    nf_run   = tuple(meta, file(".command.run", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
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
