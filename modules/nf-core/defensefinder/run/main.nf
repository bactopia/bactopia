// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'defensefinder')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::defense-finder=1.2.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
DF_VERSION     = "1.3.0"
DF_MODELS_VERSION = "1.3.0"
CASFINDER_VERSION = "3.1.0"

process DEFENSEFINDER_RUN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/defense-finder:1.3.0--pyhdfd78af_0' :
        'quay.io/biocontainers/defense-finder:1.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    each path(db)

    output:
    tuple val(meta), path("*_defense_finder_genes.tsv")  , emit: genes_tsv
    tuple val(meta), path("*_defense_finder_hmmer.tsv")  , emit: hmmer_tsv
    tuple val(meta), path("*_defense_finder_systems.tsv"), emit: systems_tsv
    tuple val(meta), path("*.prt")                       , emit: proteins
    tuple val(meta), path("*.prt.idx")                   , emit: proteins_index
    tuple val(meta), path("${prefix}.macsydata.tar.gz")  , emit: macsydata_raw, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # Extract database
    # Use custom TMPDIR to prevent FileExistsError related to writing to same tmpdir (/tmp/tmp-macsy-cache/)
    tar -xf $db
    mkdir -p defense-finder-tmp/defense-finder
    TMPDIR=defense-finder-tmp/defense-finder macsydata \\
        install \\
        --target defense-finder/ \\
        models/defense-finder-models-v${DF_MODELS_VERSION}.tar.gz

    mkdir -p defense-finder-tmp/CasFinder
    TMPDIR=defense-finder-tmp/CasFinder macsydata \\
        install \\
        --target defense-finder/ \\
        models/CasFinder-${CASFINDER_VERSION}.tar.gz

    TMPDIR=defense-finder-tmp/ defense-finder \\
        run \\
        $options.args \\
        --workers $task.cpus \\
        --models-dir defense-finder/ \\
        $fasta

    if [ "${params.df_preserveraw}" == "true" ]; then
        tar -czf ${prefix}.macsydata.tar.gz defense-finder-tmp/
        rm -rf defense-finder-tmp/ 
    fi

    # Cleanup intermediate files and unused outputs
    rm -rf models/ defense-finder/ defense-finder-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${DF_VERSION}
        defense-finder-models: ${DF_MODELS_VERSION}
        casfinder-models: ${CASFINDER_VERSION}
    END_VERSIONS
    """
}
