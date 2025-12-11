/**
 * A NATA accredited tool for reporting the presence of antimicrobial resistance genes in bacterial genomes.
 *
 * This process executes abritamr_run to perform analysis
 *
 * @status stable
 * @keywords bacteria, fasta, antibiotic resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic
 * @citation abritamr_run
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Assembled contigs in FASTA format
 *
 * @output matches   Tab-delimited file, with a row per sequence, and columns representing functional drug classes
 * @output partials  Tab-delimited file, with a row per sequence, and columns representing partial hits to functional drug classes
 * @output virulence Tab-delimited file, with a row per sequence, and columns representing AMRFinderPlus virulence gene classification
 * @output amrfinder Amrfinder
 * @output summary   Summary
 * @output logs      Optional tool execution logs
 * @output nf_logs   Nextflow execution logs
 * @output versions  Software version information (YAML format)
 */
nextflow.preview.types = true

process ABRITAMR_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta): Tuple<Map, Path>

    output:
    matches   = tuple(meta, file("${prefix}.summary_matches.txt"))
    partials  = tuple(meta, file("${prefix}.summary_partials.txt"))
    virulence = tuple(meta, file("${prefix}.summary_virulence.txt"))
    amrfinder = tuple(meta, file("${prefix}.amrfinder.out"))
    summary   = tuple(meta, file("${prefix}.abritamr.txt", optional: true))
    logs      = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs   = tuple(meta, files(".command.*"))
    versions  = tuple(meta, file("versions.yml"))

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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    abritamr run \\
        --contigs ${fasta_name} \\
        --prefix ${prefix} \\
        ${task.ext.args} \\
        --jobs ${task.cpus}

    # Rename output files to prevent name collisions
    mv ${prefix}/summary_matches.txt ./${prefix}.summary_matches.txt
    mv ${prefix}/summary_partials.txt ./${prefix}.summary_partials.txt
    mv ${prefix}/summary_virulence.txt ./${prefix}.summary_virulence.txt
    mv ${prefix}/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f ${prefix}/abritamr.txt ]; then
        # This file is not always present
        mv ${prefix}/abritamr.txt ./${prefix}.abritamr.txt
    fi

    # Cleanup
    rm -rf ${fasta_name} ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
