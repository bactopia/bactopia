nextflow.preview.types = true

process MASH_DIST {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, query) : Tuple<Map, Path>
    reference      : Path

    output:
    dist     = tuple(meta, file("*.txt"))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
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

process MERLIN_DIST {
    // Used by Merlin to extract species with matches
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, query, reads) : Tuple<Map, Path, Path>
    reference             : Path

    stage:
    stageAs 'inputs/*', query
    stageAs 'inputs/*', reads

    output:
    dist               = tuple(meta, file("*.txt"))
    escherichia        = tuple(meta, file(query, optional: true), file("escherichia.*", optional: true))
    escherichia_fq     = tuple(meta, file(reads, optional: true), file("escherichia.*", optional: true))
    escherichia_fna_fq = tuple(meta, file(query, optional: true), file(reads, optional: true), file("escherichia.*", optional: true))
    haemophilus        = tuple(meta, file(query, optional: true), file("haemophilus.*", optional: true))
    klebsiella         = tuple(meta, file(query, optional: true), file("klebsiella.*", optional: true))
    legionella         = tuple(meta, file(query, optional: true), file("legionella.*", optional: true))
    listeria           = tuple(meta, file(query, optional: true), file("listeria.*", optional: true))
    mycobacterium      = tuple(meta, file(query, optional: true), file("mycobacterium.*", optional: true))
    mycobacterium_fq   = tuple(meta, file(reads, optional: true), file("mycobacterium.*", optional: true))
    neisseria          = tuple(meta, file(query, optional: true), file("neisseria.*", optional: true))
    pseudomonas        = tuple(meta, file(query, optional: true), file("pseudomonas.*", optional: true))
    salmonella         = tuple(meta, file(query, optional: true), file("salmonella.*", optional: true))
    salmonella_fq      = tuple(meta, file(reads, optional: true), file("salmonella.*", optional: true))
    staphylococcus     = tuple(meta, file(query, optional: true), file("staphylococcus.*", optional: true))
    streptococcus      = tuple(meta, file(query, optional: true), file("streptococcus.*", optional: true))
    streptococcus_fq   = tuple(meta, file(reads, optional: true), file("streptococcus.*", optional: true))
    genus              = tuple(meta, file("*.genus", optional: true))
    logs               = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin           = tuple(meta, file(".command.begin"))
    nf_err             = tuple(meta, file(".command.err"))
    nf_log             = tuple(meta, file(".command.log"))
    nf_out             = tuple(meta, file(".command.out"))
    nf_run             = tuple(meta, file(".command.run"))
    nf_sh              = tuple(meta, file(".command.sh"))
    nf_trace           = tuple(meta, file(".command.trace"))
    versions           = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.runtype = _meta.runtype
    meta.is_compressed = _meta.is_compressed
    meta.single_end = _meta.single_end
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
        -C \\
        -p ${task.cpus} \\
        ${task.ext.args} \\
        ${reference_name} \\
        ${query} | sort -rn -k5,5 -t\$'\t' >> ${prefix}-dist.txt

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
        elif [ "${task.ext.full_merlin}" == "true" ]; then
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

    # Clean up
    rm ${reference_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
    END_VERSIONS
    """
}
