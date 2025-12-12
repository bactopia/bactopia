/**
 * Detect resistance and lineages of Mycobacterium tuberculosis genomes.
 *
 * Uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to profile *Mycobacterium tuberculosis*
 * data for drug resistance and lineage information by aligning reads to a reference genome and identifying
 * specific variants.
 *
 * @status stable
 * @keywords tuberculosis, mycobacterium, drug resistance, amr, typing, variant calling
 * @tags complexity:moderate input-type:single output-type:multiple features:compression,conditional-logic
 * @citation tbprofiler_profile
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ reads (single or paired-end)
 *
 * @output bam      Aligned BAM file
 * @output csv      Results in CSV format
 * @output json     Compressed JSON results file
 * @output txt      Results in text format
 * @output vcf      Compressed VCF file with variants
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process TBPROFILER_PROFILE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>

    output:
    bam      = tuple(meta, files("bam/*.bam"))
    csv      = tuple(meta, files("supplemental/*.csv", optional: true))
    json     = tuple(meta, files("supplemental/*.json.gz"))
    txt      = tuple(meta, files("supplemental/*.txt", optional: true))
    vcf      = tuple(meta, files("vcf/*.vcf.gz"))
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
    def input_reads = meta.single_end ? "--read1 ${reads[0]}" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    def platform = meta.runtype == "ont" ? "--platform nanopore" : "--platform illumina"
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    tb-profiler \\
        profile \\
        ${task.ext.args} \\
        ${platform} \\
        --csv \\
        --txt \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --no_trim \\
        --db_dir database/ \\
        ${input_reads}

    # Cleanup
    mv results/ supplemental/
    gzip supplemental/*.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler profile --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
