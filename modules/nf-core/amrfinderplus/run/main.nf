// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::ncbi-amrfinderplus=3.12.8"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.12.8--h283d18e_0' :
        'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0' }"

    input:
    tuple val(meta), path(genes), path(proteins), path(gff)
    each path(db)

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "*.{log,err}"                              , emit: logs, optional: true
    path ".command.*"                               , emit: nf_logs
    path "versions.yml"                             , emit: versions

    script:
    def fna_is_compressed = genes.getName().endsWith(".gz") ? true : false
    def faa_is_compressed = proteins.getName().endsWith(".gz") ? true : false
    def gff_is_compressed = gff.getName().endsWith(".gz") ? true : false
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fna_name = genes.getName().replace(".gz", "")
    faa_name = proteins.getName().replace(".gz", "")
    gff_name = gff.getName().replace(".gz", "")
    annotation_format = gff_name.endsWith(".gff") ? "prokka" : "bakta"
    """
    if [ "$fna_is_compressed" == "true" ]; then
        gzip -c -d $genes > $fna_name
    fi

    if [ "$faa_is_compressed" == "true" ]; then
        gzip -c -d $proteins > $faa_name
    fi

    if [ "$gff_is_compressed" == "true" ]; then
        gzip -c -d $gff > $gff_name
    fi

    # Extract database
    tar xzf $db

    # Full AMRFinderPlus search combining results
    amrfinder \\
        --nucleotide $fna_name \\
        --protein $faa_name \\
        --gff $gff_name \\
        --annotation_format $annotation_format \\
        $organism_param \\
        $options.args \\
        --database amrfinderplus/ \\
        --threads $task.cpus \\
        --name $prefix > ${prefix}.tsv

    # Clean up
    DB_VERSION=\$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    rm -rf amrfinderplus/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$DB_VERSION
    END_VERSIONS
    """
}
