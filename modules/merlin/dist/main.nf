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
 * @input record(meta, fna, r1, r2, se, lr)
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
 * @output record(meta, dist, results, fna, r1, r2, se, lr, escherichia, haemophilus, klebsiella, legionella, listeria, mycobacterium, neisseria, pseudomonas, salmonella, staphylococcus, streptococcus, genus, logs, nf_logs, versions)
 * - `dist`: Raw Mash distance results
 * - `fna`: Passthrough of assembled contigs
 * - `r1`: Passthrough of Illumina R1 reads
 * - `r2`: Passthrough of Illumina R2 reads
 * - `se`: Passthrough of single-end reads
 * - `lr`: Passthrough of long reads
 * - `escherichia`: Conditional marker file triggering Escherichia analysis tools
 * - `haemophilus`: Conditional marker file triggering Haemophilus analysis tools
 * - `klebsiella`: Conditional marker file triggering Klebsiella analysis tools
 * - `legionella`: Conditional marker file triggering Legionella analysis tools
 * - `listeria`: Conditional marker file triggering Listeria analysis tools
 * - `mycobacterium`: Conditional marker file triggering Mycobacterium analysis tools
 * - `neisseria`: Conditional marker file triggering Neisseria analysis tools
 * - `pseudomonas`: Conditional marker file triggering Pseudomonas analysis tools
 * - `salmonella`: Conditional marker file triggering Salmonella analysis tools
 * - `staphylococcus`: Conditional marker file triggering Staphylococcus analysis tools
 * - `streptococcus`: Conditional marker file triggering Streptococcus analysis tools
 * - `genus`: Marker file indicating the detected genus
 */
nextflow.preview.types = true

process MERLIN_DIST {
    // Used by Merlin to extract species with matches
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, fna: Path, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
    reference                                                          : Path

    stage:
    stageAs 'fna/*', fna
    stageAs 'reads/r1/*', r1
    stageAs 'reads/r2/*', r2
    stageAs 'reads/se/*', se
    stageAs 'reads/lr/*', lr

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        fna: file(fna),
        r1: file(r1 ?: 'EMPTY_R1', optional: true),
        r2: file(r2 ?: 'EMPTY_R2', optional: true),
        se: file(se ?: 'EMPTY_SE', optional: true),
        lr: file(lr ?: 'EMPTY_LR', optional: true),
        escherichia: file("escherichia.genus", optional: true),
        haemophilus: file("haemophilus.genus", optional: true),
        klebsiella: file("klebsiella.genus", optional: true),
        legionella: file("legionella.genus", optional: true),
        listeria: file("listeria.genus", optional: true),
        mycobacterium: file("mycobacterium.genus", optional: true),
        neisseria: file("neisseria.genus", optional: true),
        pseudomonas: file("pseudomonas.genus", optional: true),
        salmonella: file("salmonella.genus", optional: true),
        staphylococcus: file("staphylococcus.genus", optional: true),
        streptococcus: file("streptococcus.genus", optional: true),
        genus: files("*.genus", optional: true),
        dist: file("${prefix}-dist.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}-dist.txt")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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
