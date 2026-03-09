/**
 * Identify species to trigger genus-specific downstream analyses (Merlin).
 *
 * This is a specialized process for the [Merlin](https://bactopia.github.io/latest/bactopia-tools/merlin/)
 * workflow. It runs `mash dist` against a reference database and parses the results to detect
 * specific genera (e.g., *Salmonella*, *Staphylococcus*). Based on the detected genus, it
 * outputs data into specific channels to trigger targeted tools (e.g., finding *Salmonella* triggers Sistr).
 *
 * @status stable
 * @keywords merlin, mash, routing, logic, genus-specific, automation
 * @tags complexity:complex input-type:multiple output-type:multiple features:conditional-logic
 * @citation mash
 *
 * @input tuple(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @input reference
 * The reference Mash database to screen against
 *
 * @output dist                The raw Mash distance results
 * @output escherichia         Conditional channel triggering *Escherichia* analysis tools
 * @output escherichia_fq      Conditional channel triggering FASTQ-based *Escherichia* tools
 * @output escherichia_fna_fq  Conditional channel triggering tools requiring both FASTA and FASTQ
 * @output haemophilus         Conditional channel triggering *Haemophilus* analysis tools
 * @output klebsiella          Conditional channel triggering *Klebsiella* analysis tools
 * @output legionella          Conditional channel triggering *Legionella* analysis tools
 * @output listeria            Conditional channel triggering *Listeria* analysis tools
 * @output mycobacterium       Conditional channel triggering *Mycobacterium* analysis tools (Assembly)
 * @output mycobacterium_fq    Conditional channel triggering *Mycobacterium* analysis tools (FASTQ)
 * @output neisseria           Conditional channel triggering *Neisseria* analysis tools
 * @output pseudomonas         Conditional channel triggering *Pseudomonas* analysis tools
 * @output salmonella          Conditional channel triggering *Salmonella* analysis tools (Assembly)
 * @output salmonella_fq       Conditional channel triggering *Salmonella* analysis tools (FASTQ)
 * @output staphylococcus      Conditional channel triggering *Staphylococcus* analysis tools
 * @output streptococcus       Conditional channel triggering *Streptococcus* analysis tools (Assembly)
 * @output streptococcus_fq    Conditional channel triggering *Streptococcus* analysis tools (FASTQ)
 * @output genus               A marker file indicating the detected genus (for debugging)
 * @output logs                Optional software execution logs containing warnings/errors
 * @output nf_logs             Nextflow execution scripts and logs for debugging
 * @output versions            A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MERLIN_DIST {
    // Used by Merlin to extract species with matches
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fna, r1, r2, se, lr) : Tuple<Map, Path, Path?, Path?, Path?, Path?>
    reference                    : Path

    stage:
    stageAs 'fna/*', fna
    stageAs 'reads/r1/*', r1
    stageAs 'reads/r2/*', r2
    stageAs 'reads/se/*', se
    stageAs 'reads/lr/*', lr

    output:
    dist               = tuple(meta, file("${prefix}-dist.txt"))
    escherichia        = tuple(meta, file("escherichia.genus", optional: true), fna)
    escherichia_fq     = tuple(meta, file("escherichia.genus", optional: true), r1, r2, se, lr)
    escherichia_fna_fq = tuple(meta, file("escherichia.genus", optional: true), fna, r1, r2, se, lr)
    haemophilus        = tuple(meta, file("haemophilus.genus", optional: true), fna)
    klebsiella         = tuple(meta, file("klebsiella.genus", optional: true), fna)
    legionella         = tuple(meta, file("legionella.genus", optional: true), fna)
    listeria           = tuple(meta, file("listeria.genus", optional: true), fna)
    mycobacterium      = tuple(meta, file("mycobacterium.genus", optional: true), fna)
    mycobacterium_fq   = tuple(meta, file("mycobacterium.genus", optional: true), r1, r2, se, lr)
    neisseria          = tuple(meta, file("neisseria.genus", optional: true), fna)
    pseudomonas        = tuple(meta, file("pseudomonas.genus", optional: true), fna)
    salmonella         = tuple(meta, file("salmonella.genus", optional: true), fna)
    salmonella_fq      = tuple(meta, file("salmonella.genus", optional: true), r1, r2, se, lr)
    staphylococcus     = tuple(meta, file("staphylococcus.genus", optional: true), fna)
    streptococcus      = tuple(meta, file("streptococcus.genus", optional: true), fna)
    streptococcus_fq   = tuple(meta, file("streptococcus.genus", optional: true), r1, r2, se, lr)
    genus              = tuple(meta, files("*.genus", optional: true))
    logs               = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs            = tuple(meta, files(".command.*"))
    versions           = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.runtype = _meta.runtype
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
        ${fna} | sort -rn -k5,5 -t\$'\t' >> ${prefix}-dist.txt

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
