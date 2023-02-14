
// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'bracken')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::kraken2=2.1.2 bioconda::bracken=2.7 bioconda::fastq-scan=1.0.1 conda-forge::pigz=2.6 conda-forge::pandas=1.5.2 conda-forge::python=3.10 bioconda::krakentools=1.2 bioconda::krona=2.8.1"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process KRAKEN2_BRACKEN {
    tag "$meta.id"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, prefix:prefix, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "quay.io/bactopia/kraken2_bracken:2.2.0"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}.bracken.tsv")       , emit: tsv
    tuple val(meta), path('*classified*')                , emit: classified
    tuple val(meta), path('*unclassified*')              , emit: unclassified
    tuple val(meta), path("${prefix}.kraken2.report.txt"), emit: kraken2_report
    tuple val(meta), path("${prefix}.bracken.report.txt"), emit: bracken_report
    tuple val(meta), path("*.abundances.txt")            , emit: abundances
    tuple val(meta), path("*.krona.html")            , emit: krona
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    classified = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    def BRACKEN_VERSION = "2.7"
    def KRAKENTOOLS_VERSION = "1.2"
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
        $reads > kracken.out

    # Get read length
    if [ "${params.bracken_read_length}" == "0" ]; then
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
        READ_LENGTH="${params.bracken_read_length}"
    fi

    bracken \\
        $options.args2 \\
        -d \$KRAKEN_DB \\
        -r \$READ_LENGTH \\
        -i ${prefix}.kraken2.report.txt \\
        -w ${prefix}.bracken.report.txt \\
        -o bracken.temp

    # Sort bracken report by 'fraction_total_reads' (column 7)
    head -n 1 bracken.temp > ${prefix}.bracken.abundances.txt
    grep -v "fraction_total_reads\$" bracken.temp | sort -k 7 -rn >> ${prefix}.bracken.abundances.txt

    # Compress Kraken FASTQs
    pigz -p $task.cpus *.fastq

    # Adjust bracken to include unclassified and produce summary
    kraken-bracken-summary.py \\
        ${prefix} \\
        ${prefix}.kraken2.report.txt \\
        ${prefix}.bracken.report.txt \\
        ${prefix}.bracken.abundances.txt

    # Create a Krona report from reports
    if [ "${params.skip_krona}" == "false" ]; then
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
        rm *-krona.temp
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
