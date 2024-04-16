// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'panaroo')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::panaroo=1.4.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PANAROO_RUN {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panaroo:1.4.2--pyhdfd78af_0' :
        'quay.io/biocontainers/panaroo:1.4.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                                              , emit: results
    tuple val(meta), path("core-genome.aln.gz")                     , optional: true, emit: aln
    tuple val(meta), path("results/gene_presence_absence_roary.csv"), optional: true, emit: csv
    tuple val(meta), path("results/gene_presence_absence.csv")      , optional: true, emit: panaroo_csv
    path "*.{log,err}"                                              , optional: true, emit: logs
    path ".command.*"                                                               , emit: nf_logs
    path "versions.yml"                                                             , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gz" | xargs gunzip

    # Make FOFN of gff (Prokka) and gff3 (Bakta) files
    find gff/ -name "*.gff" -or -name "*.gff3" > gff-fofn.txt

    panaroo \\
        $options.args \\
        -t $task.cpus \\
        -o results \\
        -i gff-fofn.txt

    # Cleanup
    find . -name "*.fas" | xargs -I {} -P $task.cpus -n 1 gzip {}

    if [[ -f "results/core_gene_alignment.aln" ]]; then
        gzip results/core_gene_alignment.aln
        cp results/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' ))
    END_VERSIONS
    """
}
