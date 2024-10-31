// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'shigapass')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::shigapass=1.5.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
SHIGAPASS_VERSION = "1.5.0"

process SHIGAPASS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigapass:1.5.0--hdfd78af_0' :
        'quay.io/biocontainers/shigapass:1.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}_Flex_summary.tsv"), optional: true, emit: flex_tsv
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # ShigaPass does not accept compressed FASTA files
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Convert our genome path to a file with a path in it
    ls $fasta_name > ${fasta_name}_tmp.txt

    # Run ShigaPass
    ShigaPass.sh \\
        -l ${fasta_name}_tmp.txt \\
        $options.args \\
        -p "\$(dirname \$(which ShigaPass.sh))/../share/shigapass-${SHIGAPASS_VERSION}/db" \\
        -t $task.cpus \\
        -o ${prefix}

    # Remove the temporary file from above
    rm ${fasta_name}_tmp.txt

    # Convert to tab delimited and move to the pwd
    sed 's/;/\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv

    # Convert to tab delimited and move to the pwd
    [ ! -f ${prefix}/ShigaPass_Flex_summary.csv ] || sed 's/;/\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """
}
