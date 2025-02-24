// Import generic module functions 
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'kraken2')
options.btype = "tools"
conda_tools   = "bioconda::bactopia-teton=1.1.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-teton:1.1.1--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-teton:1.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path('*.kraken2.report.txt')              , emit: kraken2_report
    tuple val(meta), path('*.scrub.report.tsv')                , emit: scrub_report, optional: true
    tuple val(meta), path("*.${classified_naming}*.fastq.gz")  , emit: classified, optional: true
    tuple val(meta), path("*.${unclassified_naming}*.fastq.gz"), emit: unclassified, optional: true
    tuple val(meta), path("*.${classified_naming}*.fastq.gz")  , path("EMPTY_EXTRA"), emit: classified_extra, optional: true
    tuple val(meta), path("*.${unclassified_naming}*.fastq.gz"), path("EMPTY_EXTRA"), emit: unclassified_extra, optional: true
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.single_end = reads[1] == null ? true : false
    meta.is_paired = reads[1] == null ? false : true
    def paired = meta.single_end ? "" : "--paired"
    classified_naming = params.wf == "teton" || params.wf == "scrubber" || params.wf == "cleanyerreads" ? "host" : "classified"
    classified = meta.single_end ? "${prefix}.${classified_naming}.fastq"   : "${prefix}.${classified_naming}#.fastq"
    unclassified_naming = params.wf == "teton" || params.wf == "scrubber" || params.wf == "cleanyerreads" ? "scrubbed" : "unclassified"
    unclassified = meta.single_end ? "${prefix}.${unclassified_naming}.fastq" : "${prefix}.${unclassified_naming}#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
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
        $options.args \\
        $reads > /dev/null

    # If scrubbing, rename and summarize
    if [ "$unclassified_naming" == "scrubbed" ]; then
        # Rename scrubbed reads
        if [ "$meta.single_end" == "false" ]; then
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
    if [ "$is_tarball" == "true" ]; then
        rm -rf database
    fi

    if [[ "${params.keep_filtered_reads}" == "true" || "${params.wf}" == "scrubber" ]]; then
        # Compress Kraken FASTQs
        pigz -p $task.cpus *.fastq
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
