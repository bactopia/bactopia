// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'clermontyping')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::clermontyping=24.02"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CLERMONTYPING {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clermontyping:24.02--py312hdfd78af_1' :
        'quay.io/biocontainers/clermontyping:24.02--py312hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("results/${prefix}.phylogroups.txt"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    clermonTyping.sh \\
        --fasta $fasta_name \\
        --name results \\
        $options.args

    # Remove temporary files and rename outputs
    rm results/$fasta_name
    rm -rf results/db
    rm -rf results/results.R
    mv results/${fasta_name}.xml results/${prefix}.blast.xml
    mv results/${fasta_name}_mash_screen.tab results/${prefix}.mash.tsv
    mv results/results.html results/${prefix}.html

    # add column names to phylogroups file
    echo "sample<TAB>detected_genes<TAB>quadruplex_genes<TAB>quadruplex_alleles<TAB>phylogroup<TAB>mash_group" | sed 's/<TAB>/\t/g' > results/${prefix}.phylogroups.txt
    cat results/results_phylogroups.txt >> results/${prefix}.phylogroups.txt
    rm results/results_phylogroups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clermontyping: \$(echo \$(clermonTyping.sh -v 2>&1) | sed 's/^.* version : //;s/ .*\$//')
    END_VERSIONS
    """
}
