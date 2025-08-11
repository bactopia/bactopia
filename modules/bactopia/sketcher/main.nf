process SKETCHER {
    tag "${meta.id}"
    label "process_low"

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    tuple val(meta), path(fasta)
    path mash_db
    path sourmash_db

    output:
    tuple val(meta), path("${prefix}.sig")                        , emit: sig
    tuple val(meta), path("${prefix}-k*.msh")                     , emit: msh
    tuple val(meta), path("${prefix}-mash-refseq88-k21.txt")      , emit: mash
    tuple val(meta), path("${prefix}-sourmash-gtdb-rs207-k31.txt"), emit: sourmash
    path "*.{log,err}", emit: logs, optional: true
    path ".command.out", emit: nf_out
    path ".command.err", emit: nf_err
    path ".command.log", emit: nf_log
    path ".command.sh", emit: nf_sh
    path ".command.trace", emit: nf_trace
    path ".command.run", emit: nf_run, optional: true
    path ".command.begin", emit: nf_begin
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.suffix ? "${task.ext.suffix}" : "${meta.id}"
    def is_compressed = mash_db.getName().endsWith(".xz") ? true : false
    def mash_name = mash_db.getName().replace(".xz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        xz -c -d $mash_db > $mash_name
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
