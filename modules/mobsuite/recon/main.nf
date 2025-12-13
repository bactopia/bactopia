/**
 * Reconstruct and type plasmids from a bacterial genome assembly.
 *
 * Uses [MobSuite's mob_recon](https://github.com/phac-nml/mob-suite) to reconstruct plasmids
 * by grouping relevant contigs. It then uses the Mob-typer component to classify the plasmids
 * based on replicon type, incompatibility group (Inc type), and predicted mobility.
 *
 * @status stable
 * @keywords bacteria, plasmid, reconstruction, mobtyper, replicon, contigs, assembly
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation mobsuite
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output chromosome    FASTA file containing all contigs assigned to the chromosome
 * @output contig_report Assignment of each input contig to the chromosome or a specific plasmid group
 * @output plasmids      Individual FASTA files for each reconstructed plasmid group
 * @output txt           The final Mob-typer classification report (replicon, Inc type, mobility)
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MOBSUITE_RECON {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>

    output:
    chromosome    = tuple(meta, files("${prefix}-chromosome.fasta.gz"))
    contig_report = tuple(meta, files("${prefix}-contig_report.txt"))
    plasmids      = tuple(meta, files("plasmid_*.fasta.gz", optional: true))
    txt           = tuple(meta, files("${prefix}-mobtyper.txt", optional: true))
    logs          = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs       = tuple(meta, files(".command.*"))
    versions      = tuple(meta, files("versions.yml"))

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
    def is_compressed = fasta.toList()[0].getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.toList()[0].getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    mob_recon \\
        --infile ${fasta_name} \\
        ${task.ext.args} \\
        --num_threads ${task.cpus} \\
        --outdir supplemental \\
        --sample_id ${prefix}

    if [[ -f "supplemental/mobtyper_results.txt" ]]; then
        mv supplemental/mobtyper_results.txt ${prefix}-mobtyper.txt
    fi

    if [[ -f "supplemental/chromosome.fasta" ]]; then
        mv supplemental/chromosome.fasta ${prefix}-chromosome.fasta
        gzip ${prefix}-chromosome.fasta
    fi

    if [[ -f "supplemental/contig_report.txt" ]]; then
        mv supplemental/contig_report.txt ${prefix}-contig_report.txt
    fi

    # Cleanup
    gzip supplemental/*.fasta
    mv supplemental/*.fasta.gz ./
    rm -rf ${fasta_name} supplemental/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
