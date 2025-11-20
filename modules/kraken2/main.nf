nextflow.preview.types = true

process KRAKEN2 {
    tag "${prefix}"
    label 'process_high'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, Path>
    db             : Path

    output:
    kraken2_report       = tuple(meta, file('*.kraken2.report.txt'))
    scrub_report         = tuple(meta, file('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, file('*.scrub.report.tsv', optional: true))
    classified           = tuple(meta, file("*.${classified_naming}*.fastq.gz", optional: true))
    unclassified         = tuple(meta, file("*.${unclassified_naming}*.fastq.gz", optional: true))
    classified_extra     = tuple(meta, file("*.${classified_naming}*.fastq.gz", optional: true), file("EMPTY_EXTRA", optional: true))
    unclassified_extra   = tuple(meta, file("*.${unclassified_naming}*.fastq.gz", optional: true), file("EMPTY_EXTRA", optional: true))
    logs                 = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin             = tuple(meta, file(".command.begin"))
    nf_err               = tuple(meta, file(".command.err"))
    nf_log               = tuple(meta, file(".command.log"))
    nf_out               = tuple(meta, file(".command.out"))
    nf_run               = tuple(meta, file(".command.run"))
    nf_sh                = tuple(meta, file(".command.sh"))
    nf_trace             = tuple(meta, file(".command.trace"))
    versions             = tuple(meta, file("versions.yml"))

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
    meta.single_end = reads[1] == null ? true : false
    meta.is_paired = reads[1] == null ? false : true
    meta.runtype = _meta.runtype
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
        ${reads} > /dev/null

    # If scrubbing, rename and summarize
    if [ "${unclassified_naming}" == "scrubbed" ]; then
        # Rename scrubbed reads
        if [ "${meta.single_end}" == "false" ]; then
            mv ${prefix}.${unclassified_naming}_1.fastq ${prefix}_R1.scrubbed.fastq
            mv ${prefix}.${unclassified_naming}_2.fastq ${prefix}_R2.scrubbed.fastq
        fi

        # Quick stats on reads
        zcat ${reads} | fastq-scan > original.json
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
