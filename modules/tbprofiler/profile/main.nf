/**
 * Detect resistance and lineages of Mycobacterium tuberculosis genomes.
 *
 * Uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to profile *Mycobacterium tuberculosis*
 * data for drug resistance and lineage information by aligning reads to a reference genome and identifying
 * specific variants.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords tuberculosis, mycobacterium, drug resistance, amr, typing, variant calling
 * @tags complexity:moderate input-type:single output-type:multiple features:compression,conditional-logic
 * @citation tbprofiler
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
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
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>

    output:
    bam      = tuple(meta, files("bam/*.bam"))
    csv      = tuple(meta, files("supplemental/*.csv", optional: true))
    json     = tuple(meta, files("supplemental/*.json.gz"))
    txt      = tuple(meta, files("supplemental/*.txt", optional: true))
    vcf      = tuple(meta, files("vcf/*.vcf.gz"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs and platform for tb-profiler
    def input_reads = has_lr ? "--read1 ${lr}" : (meta.single_end ? "--read1 ${se}" : "--read1 ${r1} --read2 ${r2}")
    def platform = has_lr ? "--platform nanopore" : "--platform illumina"
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
