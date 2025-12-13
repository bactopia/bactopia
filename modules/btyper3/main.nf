/**
 * In silico typing and characterization of *Bacillus cereus* group genomes.
 *
 * Uses [BTyper3](https://github.com/lmc297/BTyper3) to classify *B. cereus* group isolates.
 * It determines the PanC clade, Multi-Locus Sequence Type (MLST), and screens for virulence
 * factors, crystal toxins (Bt), and antimicrobial resistance genes.
 *
 * @status stable
 * @keywords bacteria, bacillus, cereus, typing, virulence, toxin, amr, mlst
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation btyper3
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv           A tab-delimited final summary report of the typing results
 * @output supplemental  Directory containing detailed per-gene reports and raw BLAST outputs
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process BTYPER3 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>

    output:
    tsv          = tuple(meta, files("${prefix}.tsv"))
    supplemental = tuple(meta, files("supplemental/*", optional: true))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, files("versions.yml"))

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
    def is_compressed = assembly.toList()[0].getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.toList()[0].getName().replace(".gz", "")
    """
    # Btyper3 does not accept compressed files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    btyper3 \\
        ${task.ext.args} \\
        --output ./ \\
        --input ${assembly_name}

    mv btyper3_final_results/ supplemental/
    mv supplemental/${prefix}_final_results.txt ./${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name} ${assembly_name}.njs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
