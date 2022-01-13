// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'mashdist')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::mash=2.3"
conda_env   = file("${params.condadir}/mash").exists() ? "${params.condadir}/mash" : conda_tools

process MASH_DIST {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(query)
    path reference

    output:
    tuple val(meta), path("*.txt")          , emit: dist
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    echo "reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes" | sed 's/<TAB>/\t/g' > ${prefix}-dist.txt
    mash \\
        dist \\
        -p $task.cpus \\
        $options.args \\
        $reference \\
        $query >> ${prefix}-dist.txt 2> mash.stderr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}

process MERLIN_DIST {
    // Used by Merlin to extract species with matches
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(query), path(reads)
    path reference

    output:
    tuple val(meta), path("*.txt")          , emit: dist
    tuple val(meta), path(query), path("escherichia.*")   , emit: escherichia, optional: true
    tuple val(meta), path(query), path("haemophilus.*")   , emit: haemophilus, optional: true
    tuple val(meta), path(query), path("klebsiella.*")    , emit: klebsiella, optional: true
    tuple val(meta), path(query), path("listeria.*")      , emit: listeria, optional: true
    tuple val(meta), path(query), path("mycobacterium.*") , emit: mycobacterium, optional: true
    tuple val(meta), path(reads), path("mycobacterium.*") , emit: mycobacterium_fq, optional: true
    tuple val(meta), path(query), path("neisseria.*")     , emit: neisseria, optional: true
    tuple val(meta), path(query), path("salmonella.*")    , emit: salmonella, optional: true
    tuple val(meta), path(query), path("staphylococcus.*"), emit: staphylococcus, optional: true
    tuple val(meta), path(query), path("streptococcus.*") , emit: streptococcus, optional: true
    path "*.genus"                          , emit: genus, optional: true
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes" | sed 's/<TAB>/\t/g' > ${prefix}-dist.txt
    mash \\
        dist \\
        -C \\
        -p $task.cpus \\
        $options.args \\
        $reference \\
        $query | sort -rn -k5,5 -t\$'\t' >> ${prefix}-dist.txt 2> mash.stderr.txt

    # Extract genus with hits
    declare -a GENUS=(
        "escherichia" "haemophilus" "klebsiella" "listeria" "mycobacterium" "neisseria" "salmonella" "shigella" "staphylococcus" "streptococcus"
    )
    for i in "\${GENUS[@]}"; do
        if grep -q -i "\${i}" ${prefix}-dist.txt; then
            if [ "\${i}" == "shigella" ]; then
                touch escherichia.genus
            else
                touch \${i}.genus
            fi
        elif [ "${params.full_merlin}" == "true" ]; then
            if [ "\${i}" == "shigella" ]; then
                touch escherichia.genus
            else
                if [ "\${i}" != "listeria" ]; then
                    # lissero fails on non-Listeria samples
                    touch \${i}.genus
                fi
            fi
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
