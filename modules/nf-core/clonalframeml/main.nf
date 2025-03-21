// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'clonalframeml')
options.btype = "comparative"
conda_tools   = "bioconda::clonalframeml=1.12 bioconda::maskrc-svg=0.5"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CLONALFRAMEML {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f5c68f1508671d5744655da9b0e8b609098f4138:7e089189af7822a6a18245830639dbfe11a4c277-0' :
        'quay.io/biocontainers/mulled-v2-f5c68f1508671d5744655da9b0e8b609098f4138:7e089189af7822a6a18245830639dbfe11a4c277-0' }"

    input:
    tuple val(meta), path(msa), path(newick)

    output:
    tuple val(meta), path("*.emsim.txt")                     , emit: emsim, optional: true
    tuple val(meta), path("*.em.txt")                        , emit: em
    tuple val(meta), path("*.importation_status.txt")        , emit: status
    tuple val(meta), path("*.labelled_tree.newick")          , emit: newick
    tuple val(meta), path("*.ML_sequence.fasta.gz")          , emit: fasta
    tuple val(meta), path("*.position_cross_reference.txt.gz"), emit: pos_ref
    tuple val(meta), path("*.masked.aln.gz")                 , emit: masked_aln
    path "*.{log,err}"                                       , emit: logs, optional: true
    path ".command.*"                                        , emit: nf_logs
    path "versions.yml"                                      , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    ClonalFrameML \\
        $newick \\
        $msa_name \\
        $prefix \\
        $options.args

    maskrc-svg.py $prefix --aln ${msa_name} --symbol '-' --out ${prefix}.masked.aln
    gzip ${prefix}.masked.aln

    # Cleanup
    rm ${msa_name}
    gzip *.ML_sequence.fasta
    gzip *.position_cross_reference.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clonalframeml: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
        maskrc-svg: \$( echo \$(maskrc-svg.py --version 2>&1) | sed 's/^.*maskrc-svg.py //' )
    END_VERSIONS
    """
}
