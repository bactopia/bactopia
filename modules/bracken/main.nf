process BRACKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}.bracken.tsv")                    , emit: tsv
    tuple val(meta), path('*classified*')                             , emit: classified, optional: true
    tuple val(meta), path('*unclassified*')                           , emit: unclassified, optional: true
    tuple val(meta), path("${prefix}.kraken2.report.txt")             , emit: kraken2_report
    tuple val(meta), path("${prefix}.kraken2.output.txt")             , emit: kraken2_output, optional: true
    tuple val(meta), path("${prefix}.bracken.report.txt")             , emit: bracken_report
    tuple val(meta), path("*.krona.html")                             , emit: krona, optional: true
    tuple val(meta), path("${prefix}.bracken.abundances.txt")         , emit: abundances
    tuple val(meta), path("${prefix}.bracken.classification.txt")     , emit: classification
    tuple val(meta), path("${prefix}.bracken.adjusted.abundances.txt"), emit: adjusted_abundances
    tuple val(meta), path("*.{log,err}" )  , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def paired = meta.single_end ? "" : "--paired"
    classified = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    def BRACKEN_VERSION = "2.7"
    def KRAKENTOOLS_VERSION = "1.2"
    meta.teton_reads = reads.join(",")
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        KRAKEN_DB=\$(find database/ -name "hash.k2d" | sed 's=hash.k2d==')
    else
        KRAKEN_DB=\$(find $db/ -name "hash.k2d" | sed 's=hash.k2d==')
    fi

    kraken2 \\
        --db \$KRAKEN_DB \\
        --threads $task.cpus \\
        --unclassified-out $unclassified \\
        --classified-out $classified \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $paired \\
        $task.ext.args \\
        $reads > ${prefix}.kraken2.output.txt

    # Get read length
    if [ "${task.ext.bracken_read_length}" == "0" ]; then
        OBS_READ_LENGTH=\$(zcat ${reads[0]} | fastq-scan -q | jq -r '.qc_stats.read_median')
        echo \$OBS_READ_LENGTH
        # Pre-built Bracken databases come with 50,75,100,150,200,250,300, split the difference
        if [ "\$OBS_READ_LENGTH" -gt 275 ]; then
            READ_LENGTH="300"
        elif [ "\$OBS_READ_LENGTH" -gt 225 ]; then
            READ_LENGTH="250"
        elif [ "\$OBS_READ_LENGTH" -gt 175 ]; then
            READ_LENGTH="200"
        elif [ "\$OBS_READ_LENGTH" -gt 125 ]; then
            READ_LENGTH="150"
        elif [ "\$OBS_READ_LENGTH" -gt 85 ]; then
            READ_LENGTH="100"
        elif [ "\$OBS_READ_LENGTH" -gt 65 ]; then
            READ_LENGTH="75"
        else
            READ_LENGTH="50"
        fi
    else
        # use user defined read length
        READ_LENGTH="${task.ext.bracken_read_length}"
    fi

    bracken \\
        $task.ext.args2 \\
        -d \$KRAKEN_DB \\
        -r \$READ_LENGTH \\
        -i ${prefix}.kraken2.report.txt \\
        -w ${prefix}.bracken.report.txt \\
        -o bracken.temp

    # Sort bracken report by 'fraction_total_reads' (column 7)
    head -n 1 bracken.temp > ${prefix}.bracken.abundances.txt
    grep -v "fraction_total_reads\$" bracken.temp | sort -k 7 -rn -t \$'\t' >> ${prefix}.bracken.abundances.txt

    # Adjust bracken to include unclassified and produce summary
    kraken-bracken-summary.py \\
        ${prefix} \\
        ${prefix}.kraken2.report.txt \\
        ${prefix}.bracken.report.txt \\
        ${prefix}.bracken.abundances.txt \\
        --max_secondary_percent ${task.ext.bracken_max_secondary_percent}

    # Create a Krona report from reports
    if [ "${task.ext.skip_krona}" == "false" ]; then
        # Kraken2
        kreport2krona.py \\
            --report ${prefix}.kraken2.report.txt \\
            --output kraken2-krona.temp
        ktImportText -o ${prefix}.kraken2.krona.html kraken2-krona.temp

        # Bracken
        kreport2krona.py \\
            --report ${prefix}.bracken.report.txt \\
            --output bracken-krona.temp
        ktImportText -o ${prefix}.bracken.krona.html bracken-krona.temp
    fi

    # Clean up large files produced by Kraken2/Bracken
    rm *.temp
    if [ "${task.ext.kraken2_keep_raw_output}" == "false" ]; then
        # Remove kraken2 STDOUT output file
        rm ${prefix}.kraken2.output.txt
    fi

    if [ "${task.ext.kraken2_keep_filtered_reads}" == "true" ]; then
        # Compress Kraken FASTQs
        pigz -p $task.cpus *.fastq
    else
        # Remove filtered FASTQs
        rm *.fastq
    fi

    if [ "$is_tarball" == "true" ]; then
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${BRACKEN_VERSION}
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        jq: \$(echo \$(jq --version 2>&1) | sed 's/jq-//')
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        krakentools: ${KRAKENTOOLS_VERSION}
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
