/**
 * Taxonomic classification and host filtering of sequence reads.
 *
 * Uses [Kraken2](https://github.com/DerrickWood/kraken2) to assign taxonomic labels to short
 * DNA reads by examining exact k-mer matches against a large reference database. It uses the
 * Lowest Common Ancestor (LCA) algorithm to provide high-precision classification, making it
 * ideal for metagenomics or removing host contamination (scrubbing).
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomy, classification, contamination, scrubbing, k-mer, lca
 * @tags complexity:complex input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation kraken2
 *
 * @note Database Required
 * Requires a standard Kraken2 database (directory or tarball). Memory usage depends on database size (Standard ~50GB).
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) - not typically used by Kraken2
 *
 * @input db
 * Kraken2 database (Directory or compressed tarball)
 *
 * @output kraken2_report    Standard Kraken2 report containing taxonomic abundance counts
 * @output scrub_report      Summary report of reads removed during host scrubbing (optional)
 * @output scrub_special_report Duplicate scrub report with modified metadata (Internal use)
 * @output classified        Reads assigned to a taxon in the database (FASTQ)
 * @output unclassified      Reads NOT assigned to any taxon (FASTQ)
 * @output classified_extra  Duplicate classified channel with placeholder for pipeline routing
 * @output unclassified_extra Duplicate unclassified channel with placeholder for pipeline routing
 * @output logs              Optional software execution logs containing warnings/errors
 * @output nf_logs           Nextflow execution scripts and logs for debugging
 * @output versions          A YAML formatted file with software versions
 */
nextflow.preview.types = true

process KRAKEN2 {
    tag "${prefix}"
    label 'process_high'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>
    db                      : Path

    output:
    kraken2_report       = tuple(meta, files('*.kraken2.report.txt'))
    scrub_report         = tuple(meta, files('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, files('*.scrub.report.tsv', optional: true))
    classified           = tuple(meta, files("*.${classified_naming}*.fastq.gz", optional: true))
    unclassified         = tuple(meta, files("*.${unclassified_naming}*.fastq.gz", optional: true))
    classified_extra     = tuple(meta, files("*.${classified_naming}*.fastq.gz", optional: true), files("EMPTY_EXTRA", optional: true))
    unclassified_extra   = tuple(meta, files("*.${unclassified_naming}*.fastq.gz", optional: true), files("EMPTY_EXTRA", optional: true))
    logs                 = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs              = tuple(meta, files(".command.*"))
    versions             = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    output_folder = task.ext.wf == "scrubber" || task.ext.wf == "teton" ? "scrubber" : "${task.ext.process_name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope

    if (task.ext.wf == "teton") {
        meta.output_dir = "${prefix}/teton/tools/${output_folder}"
        meta.logs_dir = "${prefix}/teton/tools/${output_folder}/logs/${task.ext.logs_subdir}"
    }
    else {
        meta.output_dir = "${prefix}/tools/${output_folder}"
        meta.logs_dir = "${prefix}/tools/${output_folder}/logs/${task.ext.logs_subdir}"
    }
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2
    meta.runtype = _meta.runtype

    // Build read inputs for kraken2
    read_inputs = meta.single_end ? "${se}" : "${r1} ${r2}"

    special_meta = [:]
    special_meta.id = prefix
    def paired = meta.single_end ? "" : "--paired"
    classified_naming = task.ext.wf != "kraken2" ? "host" : "classified"
    classified = meta.single_end ? "${prefix}.${classified_naming}.fastq" : "${prefix}.${classified_naming}#.fastq"
    unclassified_naming = task.ext.wf != "kraken2" ? "scrubbed" : "unclassified"
    unclassified = meta.single_end ? "${prefix}.${unclassified_naming}.fastq" : "${prefix}.${unclassified_naming}#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        KRAKEN_DB=\$(find database/ -name "hash.k2d" | sed 's=hash.k2d==')
    else
        KRAKEN_DB=\$(find ${db}/ -name "hash.k2d" | sed 's=hash.k2d==')
    fi

    kraken2 \\
        --db \$KRAKEN_DB \\
        --threads ${task.cpus} \\
        --unclassified-out ${unclassified} \\
        --classified-out ${classified} \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        ${paired} \\
        ${task.ext.args} \\
        ${read_inputs} > /dev/null

    # If scrubbing, rename and summarize
    if [ "${unclassified_naming}" == "scrubbed" ]; then
        # Rename scrubbed reads
        if [ "${meta.single_end}" == "false" ]; then
            mv ${prefix}.${unclassified_naming}_1.fastq ${prefix}_R1.scrubbed.fastq
            mv ${prefix}.${unclassified_naming}_2.fastq ${prefix}_R2.scrubbed.fastq
        fi

        # Quick stats on reads
        zcat ${read_inputs} | fastq-scan > original.json
        cat *.scrubbed.fastq | fastq-scan > scrubbed.json
        scrubber-summary.py ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Remove host reads and temp json files
        rm ${prefix}.host*.fastq original.json scrubbed.json
    fi

    # Clean up database and large files produced by Kraken2
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database
    fi

    if [[ "${task.ext.keep_filtered_reads}" == "true" || "${task.ext.wf}" == "scrubber" || "${task.ext.wf}" == "teton" ]]; then
        # Compress Kraken FASTQs
        pigz -p ${task.cpus} *.fastq
    else
        # Remove filtered FASTQs
        rm *.fastq
    fi

    # Used for clean-yer-reads
    touch EMPTY_EXTRA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
