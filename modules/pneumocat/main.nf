/**
 * Capsular typing of Streptococcus pneumoniae from Illumina reads.
 *
 * Uses [PneumoCaT](https://github.com/ukhsa-collaboration/PneumoCaT) (Pneumococcal Capsular Typing)
 * to assign capsular types to *Streptococcus pneumoniae* using a two-step approach: first matching
 * reads to a global database, then using a mapped-based approach for specific serogroup differentiation.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2) where each read slot is Path
 *
 * @status stable
 * @keywords pneumocat, streptococcus pneumoniae, capsular typing, serotyping
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation pneumocat
 *
 * @note
 * Negative results will cause non-0 exit codes from PneumoCaT
 *
 * @input record(meta, r1, r2)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 *
 * @output record(meta, xml?, txt?, results, logs, nf_logs, versions)
 * - `xml?`: The PneumoCaT result files in XML format
 * - `txt?`: A file containing the coverage information across the genes
 */
nextflow.preview.types = true

process PNEUMOCAT {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, r1: Path, r2: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        xml: file("${prefix}.xml", optional: true),
        txt: file("${prefix}.coverage_summary.txt", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.xml", optional: true),
            files("${prefix}.coverage_summary.txt", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.2.1'
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        --threads ${task.cpus} \\
        --output_dir ./

    # Cleanup
    rm -rf *.bam *.bai ComponentComplete.txt

    # PneumoCAT uses first match in a glob, so moves between R1 and R2
    if [ -f ${prefix}_R1.results.xml ]; then
        mv ${prefix}_R1.results.xml ${prefix}.xml
    else
        mv ${prefix}_R2.results.xml ${prefix}.xml
    fi
    mv logs/* ./
    if [ -f pneumo_capsular_typing.stderr ]; then
        mv pneumo_capsular_typing.stderr pneumo_capsular_typing.stderr.log
    fi
    if [ -f pneumo_capsular_typing.stdout ]; then
        mv pneumo_capsular_typing.stdout pneumo_capsular_typing.stdout.log
    fi
    mv coverage_summary.txt ${prefix}.coverage_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: ${VERSION}
    END_VERSIONS
    """
}
