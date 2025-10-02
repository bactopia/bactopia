process SKETCHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)
    path mash_db
    path sourmash_db

    output:
    tuple val(meta), path("${prefix}.sig")                        , emit: sig
    tuple val(meta), path("${prefix}-k*.msh")                     , emit: msh
    tuple val(meta), path("${prefix}-mash-refseq88-k21.txt")      , emit: mash
    tuple val(meta), path("${prefix}-sourmash-gtdb-rs207-k31.txt"), emit: sourmash
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path(".command.run")  , emit: nf_run, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
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
