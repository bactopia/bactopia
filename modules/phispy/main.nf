/**
 * Predict prophage regions integrated into bacterial genomes.
 *
 * Uses [PhiSpy](https://github.com/linsalrob/PhiSpy) to identify integrated bacteriophage
 * (prophage) regions in a fully annotated bacterial genome. The prediction relies on scoring
 * features like strand-switch, AT-skew, unique phage-like proteins, and short coding regions.
 *
 * @status stable
 * @keywords genomics, virus, phage, prophage, bacteriophage, identification, annotation
 * @tags complexity:complex input-type:single output-type:multiple
 * @citation phispy
 *
 * @input tuple(meta, gbk)
 * - `meta`: Groovy Map containing sample information
 * - `gbk`: Annotated genome file in GenBank (*.gbk or *.gbff) format
 *
 * @output tsv            Coordinates (start/end) of each predicted prophage region in the genome
 * @output information    Detailed information and confidence scores for each prophage
 * @output bacteria_fasta FASTA sequence of the host genome with prophages excised
 * @output bacteria_gbk   GenBank file of the host genome with prophages excised
 * @output phage_fasta    FASTA sequences of the excised prophage elements
 * @output phage_gbk      GenBank file of the excised prophage elements
 * @output prophage_gff   Prophage predictions in GFF3 format
 * @output prophage_tbl   Prophage predictions in table format (GenBank annotation style)
 * @output prophage_tsv   Prophage predictions in TSV format
 * @output logs           Optional software execution logs containing warnings/errors
 * @output nf_logs        Nextflow execution scripts and logs for debugging
 * @output versions       A YAML formatted file with software versions
 */
nextflow.preview.types = true

process PHISPY {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, gbk) : Tuple<Map, Path>

    output:
    tsv            = tuple(meta, file("${prefix}.tsv"))
    information    = tuple(meta, file("supplemental/${prefix}_prophage_information.tsv", optional: true))
    bacteria_fasta = tuple(meta, file("supplemental/${prefix}_bacteria.fasta", optional: true))
    bacteria_gbk   = tuple(meta, file("supplemental/${prefix}_bacteria.gbk", optional: true))
    phage_fasta    = tuple(meta, file("supplemental/${prefix}_phage.fasta", optional: true))
    phage_gbk      = tuple(meta, file("supplemental/${prefix}_phage.gbk", optional: true))
    prophage_gff   = tuple(meta, file("supplemental/${prefix}_prophage.gff3", optional: true))
    prophage_tbl   = tuple(meta, file("supplemental/${prefix}_prophage.tbl", optional: true))
    prophage_tsv   = tuple(meta, file("supplemental/${prefix}_prophage.tsv", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, file("versions.yml"))

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
