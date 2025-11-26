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
