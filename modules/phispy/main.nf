/**
 * Predict prophage regions integrated into bacterial genomes.
 *
 * Uses [PhiSpy](https://github.com/linsalrob/PhiSpy) to identify integrated bacteriophage
 * (prophage) regions in a fully annotated bacterial genome. The prediction relies on scoring
 * features like strand-switch, AT-skew, unique phage-like proteins, and short coding regions.
 *
 * @status stable
 * @keywords genomics, virus, phage, prophage, bacteriophage, identification, annotation
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic
 * @citation phispy
 *
 * @input tuple(meta, gbk)
 * - `meta`: Groovy Map containing sample information
 * - `gbk`: Annotated genome file in GenBank (*.gbk or *.gbff) format
 *
 * @output record(meta, tsv, supplemental, results, logs, nf_logs, versions)
 * - `tsv`: Coordinates (start/end) of each predicted prophage region in the genome
 * - `supplemental`: Directory containing detailed prophage information, sequences, and annotations
 */
nextflow.preview.types = true

process PHISPY {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, gbk: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("supplemental/*", optional: true)
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
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    mkdir supplemental/
    PhiSpy.py \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        -p ${prefix} \\
        -o supplemental \\
        ${gbk}

    mv supplemental/${prefix}_prophage_coordinates.tsv ${prefix}.tsv
    mv supplemental/${prefix}_phispy.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """
}
