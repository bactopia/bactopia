// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'amrfinderplus')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::ncbi-amrfinderplus=3.10.21"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.21--h17dc2d4_0':
        'quay.io/biocontainers/ncbi-amrfinderplus:3.10.21--h17dc2d4_0' }"

    input:
    tuple val(meta), path(genes), path(proteins)
    each path(db)

    output:
    tuple val(meta), path("${prefix}-genes.tsv")                     , emit: gene_report
    tuple val(meta), path("${prefix}-proteins.tsv")                  , emit: protein_report
    tuple val(meta), path("${prefix}-{genes,proteins}-mutations.tsv"), emit: mutation_reports, optional: true
    path "*.{log,err}"                                               , emit: logs, optional: true
    path ".command.*"                                                , emit: nf_logs
    path "versions.yml"                                              , emit: versions

    script:
    def fna_is_compressed = genes.getName().endsWith(".gz") ? true : false
    def faa_is_compressed = proteins.getName().endsWith(".gz") ? true : false
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    fna_organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-genes-mutations.tsv" : ""
    faa_organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-proteins-mutations.tsv" : ""
    fna_name = genes.getName().replace(".gz", "")
    faa_name = proteins.getName().replace(".gz", "")
    """
    if [ "$fna_is_compressed" == "true" ]; then
        gzip -c -d $genes > $fna_name
    fi

    if [ "$faa_is_compressed" == "true" ]; then
        gzip -c -d $proteins > $faa_name
    fi

    mkdir amrfinderdb
    tar xzvf $db -C amrfinderdb

    # Gene
    amrfinder \\
       -n $fna_name \\
        $fna_organism_param \\
        $options.args \\
        --plus \\
        --database amrfinderdb/ \\
        --threads $task.cpus \\
        --name $prefix > ${prefix}-genes.tsv

    # Protein
    amrfinder \\
       -p $faa_name \\
        $faa_organism_param \\
        $options.args \\
        --plus \\
        --database amrfinderdb/ \\
        --threads $task.cpus \\
        --name $prefix > ${prefix}-proteins.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
