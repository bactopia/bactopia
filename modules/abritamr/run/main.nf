/**
 * Detect antimicrobial resistance and virulence genes.
 *
 * Uses [abriTAMR](https://github.com/MDU-PHL/abritamr), a NATA (National Association of 
 * Testing Authorities) accredited pipeline, to report the presence of reportable AMR 
 * genes. It acts as a wrapper for AMRFinderPlus, formatted for clinical reporting standards 
 * used in Victoria, Australia.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, antimicrobial resistance, nata, amrfinderplus
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output,compression,conditional-logic
 * @citation abritamr
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output summary   A tab-delimited summary of detected resistance genes
 * @output matches   A tab-delimited file of sequences vs functional drug classes
 * @output partials  A tab-delimited file of partial hits to functional drug classes
 * @output virulence A tab-delimited file of AMRFinderPlus virulence gene classifications
 * @output amrfinder Raw output from AMRFinderPlus (per sequence)
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process ABRITAMR_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly): Tuple<Map, Path>

    output:
    summary   = tuple(meta, file("${prefix}.abritamr.tsv", optional: true))
    matches   = tuple(meta, file("${prefix}.summary_matches.tsv"))
    partials  = tuple(meta, file("${prefix}.summary_partials.tsv"))
    virulence = tuple(meta, file("${prefix}.summary_virulence.tsv"))
    amrfinder = tuple(meta, file("${prefix}.amrfinder.out"))
    logs      = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs   = tuple(meta, files(".command.*"))
    versions  = tuple(meta, files("versions.yml"))

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

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    abritamr run \\
        --contigs ${assembly_name} \\
        --prefix ${prefix} \\
        ${task.ext.args} \\
        --jobs ${task.cpus}

    # Rename output files to prevent name collisions
    mv ${prefix}/summary_matches.txt ./${prefix}.summary_matches.tsv
    mv ${prefix}/summary_partials.txt ./${prefix}.summary_partials.tsv
    mv ${prefix}/summary_virulence.txt ./${prefix}.summary_virulence.tsv
    mv ${prefix}/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f ${prefix}/abritamr.txt ]; then
        # This file is not always present
        mv ${prefix}/abritamr.txt ./${prefix}.abritamr.tsv
    fi

    # Cleanup
    rm -rf ${assembly_name} ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
