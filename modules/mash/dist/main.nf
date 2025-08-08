process MASH_DIST {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(query)
    path reference

    output:
    tuple val(meta), path("*.txt")          , emit: dist
    tuple val(meta), path("*.{log,err}")    , emit: logs, optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    prefix = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = reference.getName().endsWith(".xz") ? true : false
    def reference_name = reference.getName().replace(".xz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        xz -c -d $reference > $reference_name
    fi

    echo "reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes" | sed 's/<TAB>/\t/g' > ${prefix}-dist.txt
    mash \\
        dist \\
        -p $task.cpus \\
        $task.ext.args \\
        $reference_name \\
        $query >> ${prefix}-dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
    END_VERSIONS
    """
}

process MERLIN_DIST {
    // Used by Merlin to extract species with matches
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(query), path(reads)
    path reference

    output:
    tuple val(meta), path("*.txt"), emit: dist
    tuple val(meta), path(query), path("escherichia.*")   , emit: escherichia, optional: true
    tuple val(meta), path(reads), path("escherichia.*")   , emit: escherichia_fq, optional: true
    tuple val(meta), path(query), path(reads), path("escherichia.*"), emit: escherichia_fna_fq, optional: true
    tuple val(meta), path(query), path("haemophilus.*")   , emit: haemophilus, optional: true
    tuple val(meta), path(query), path("klebsiella.*")    , emit: klebsiella, optional: true
    tuple val(meta), path(query), path("legionella.*")    , emit: legionella, optional: true
    tuple val(meta), path(query), path("listeria.*")      , emit: listeria, optional: true
    tuple val(meta), path(query), path("mycobacterium.*") , emit: mycobacterium, optional: true
    tuple val(meta), path(reads), path("mycobacterium.*") , emit: mycobacterium_fq, optional: true
    tuple val(meta), path(query), path("neisseria.*")     , emit: neisseria, optional: true
    tuple val(meta), path(query), path("pseudomonas.*")   , emit: pseudomonas, optional: true
    tuple val(meta), path(query), path("salmonella.*")    , emit: salmonella, optional: true
    tuple val(meta), path(reads), path("salmonella.*")    , emit: salmonella_fq, optional: true
    tuple val(meta), path(query), path("staphylococcus.*"), emit: staphylococcus, optional: true
    tuple val(meta), path(query), path("streptococcus.*") , emit: streptococcus, optional: true
    tuple val(meta), path(reads), path("streptococcus.*") , emit: streptococcus_fq, optional: true
    path "*.genus"                                        , emit: genus, optional: true
    tuple val(meta), path("*.{log,err}")                 , emit: logs, optional: true
    tuple val(meta), path(".command.begin")              , emit: nf_begin
    tuple val(meta), path(".command.err")                , emit: nf_err
    tuple val(meta), path(".command.log")                , emit: nf_log
    tuple val(meta), path(".command.out")                , emit: nf_out
    tuple val(meta), path(".command.run")                , emit: nf_run
    tuple val(meta), path(".command.sh")                 , emit: nf_sh
    tuple val(meta), path(".command.trace")              , emit: nf_trace
    tuple val(meta), path("versions.yml")                , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = reference.getName().endsWith(".xz") ? true : false
    def reference_name = reference.getName().replace(".xz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        xz -c -d $reference > $reference_name
    fi

    echo "reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes" | sed 's/<TAB>/\t/g' > ${prefix}-dist.txt
    mash \\
        dist \\
        -C \\
        -p $task.cpus \\
        $task.ext.args \\
        $reference_name \\
        $query | sort -rn -k5,5 -t\$'\t' >> ${prefix}-dist.txt

    # Extract genus with hits
    declare -a GENUS=(
        "escherichia" "haemophilus" "glaesserella" "klebsiella" "legionella" "listeria" "mycobacterium" "neisseria" "pseudomonas" "salmonella" "shigella" "staphylococcus" "streptococcus"
    )
    for i in "\${GENUS[@]}"; do
        if grep -q -i "\${i}" ${prefix}-dist.txt; then
            if [ "\${i}" == "shigella" ]; then
                touch escherichia.genus
            elif [ "\${i}" == "glaesserella" ]; then
                touch haemophilus.genus
            elif [ "\${i}" == "streptococcus" ]; then
                touch streptococcus.genus
            else
                touch \${i}.genus
            fi
        elif [ "${params.full_merlin}" == "true" ]; then
            if [ "\${i}" == "shigella" ]; then
                touch escherichia.genus
            elif [ "\${i}" == "glaesserella" ]; then
                touch haemophilus.genus
            else
                if [ "\${i}" != "listeria" ]; then
                    # lissero fails on non-Listeria samples
                    touch \${i}.genus
                fi
            fi
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
    END_VERSIONS
    """
}
