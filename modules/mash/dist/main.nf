/**
 * Calculate Mash distances between reference and query seqeunces.
 *
 * This process executes mash_dist to perform analysis
 *
 * @status stable
 * @keywords mash/dist
 * @tags complexity:moderate input-type:multiple output-type:single features:compression, conditional-logic
 * @citation mash_dist
 *
 * @input tuple(meta, query)
 * - `meta`: Groovy Map containing sample information
 * - `query`: FASTA, FASTQ or Mash sketch
 *
 * @input reference
 * FASTA, FASTQ or Mash sketch
 *
 * @output dist     The results from mash dist
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
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
    dist     = tuple(meta, files("*.txt"))
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
    (_meta, query, reads) : Tuple<Map, Set<Path>, Set<Path>>
    reference             : Path

    stage:
    stageAs 'inputs/*', query
    stageAs 'inputs/*', reads

    output:
    dist               = tuple(meta, files("*.txt"))
    escherichia        = tuple(meta, query, file("escherichia.genus", optional: true))
    escherichia_fq     = tuple(meta, reads, file("escherichia.genus", optional: true))
    escherichia_fna_fq = tuple(meta, query, reads, file("escherichia.genus", optional: true))
    haemophilus        = tuple(meta, query, file("haemophilus.genus", optional: true))
    klebsiella         = tuple(meta, query, file("klebsiella.genus", optional: true))
    legionella         = tuple(meta, query, file("legionella.genus", optional: true))
    listeria           = tuple(meta, query, file("listeria.genus", optional: true))
    mycobacterium      = tuple(meta, query, file("mycobacterium.genus", optional: true))
    mycobacterium_fq   = tuple(meta, reads, file("mycobacterium.genus", optional: true))
    neisseria          = tuple(meta, query, file("neisseria.genus", optional: true))
    pseudomonas        = tuple(meta, query, file("pseudomonas.genus", optional: true))
    salmonella         = tuple(meta, query, file("salmonella.genus", optional: true))
    salmonella_fq      = tuple(meta, reads, file("salmonella.genus", optional: true))
    staphylococcus     = tuple(meta, query, file("staphylococcus.genus", optional: true))
    streptococcus      = tuple(meta, query, file("streptococcus.genus", optional: true))
    streptococcus_fq   = tuple(meta, reads, file("streptococcus.genus", optional: true))
    genus              = tuple(meta, files("*.genus", optional: true))
    logs               = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs            = tuple(meta, files(".command.*"))
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
