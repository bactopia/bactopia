// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'hicap')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::hicap=1.0.4"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process HICAP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicap:1.0.4--pyhdfd78af_1' :
        'quay.io/biocontainers/hicap:1.0.4--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path database_dir
    path model_fp

    output:
    tuple val(meta), path("*.gbk"), emit: gbk, optional: true
    tuple val(meta), path("*.svg"), emit: svg, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    hicap \\
        --query_fp $fasta_name \\
        $database_args \\
        $model_args \\
        $options.args \\
        --threads $task.cpus \\
        --debug \\
        -o ./

    if [ ! -f ${meta.id}.tsv ]; then
        echo "isolate<TAB>predicted_serotype<TAB>attributes<TAB>genes_identified<TAB>locus_location<TAB>region_I_genes<TAB>region_II_genes<TAB>region_III_genes<TAB>IS1016_hits" | sed 's/<TAB>/\t/g' > ${meta.id}.tsv
        echo "${meta.id}<TAB>cap_not_found<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-" | sed 's/<TAB>/\t/g' >> ${meta.id}.tsv
    else
        sed -i 's/#isolate/isolate/' ${meta.id}.tsv
    fi

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
