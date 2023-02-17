// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'fastani')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::fastani=1.33"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process FASTANI {
    tag "${reference_name}"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.33--h0fdf51a_1' :
        'quay.io/biocontainers/fastani:1.33--h0fdf51a_1' }"


    input:
    tuple val(meta), path(query, stageAs: 'query-tmp/*')
    each path(reference)

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    reference_fasta = reference.getName().replace(".gz", "")
    reference_name = reference_fasta.replace(".fna", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $reference_fasta
    fi

    mkdir query
    cp -L query-tmp/* query/
    find query/ -name "*.gz" | xargs gunzip
    find query/ -name "*" -type f > query-list.txt

    fastANI \\
        --ql query-list.txt \\
        -r $reference_fasta \\
        -o fastani-result.tmp

    echo "query<TAB>reference<TAB>ani<TAB>mapped_fragments<TAB>total_fragments" | sed 's/<TAB>/\t/g' > ${reference_name}.tsv
    sed 's=^query/==' fastani-result.tmp>> ${reference_name}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
