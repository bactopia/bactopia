// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'phyloflash')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::phyloflash=3.4.2"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PHYLOFLASH  {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phyloflash:3.4.2--hdfd78af_0' :
        'quay.io/biocontainers/phyloflash:3.4.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}/*")     , emit: results
    file "${sample}/${sample}.toalign.fasta"  , emit: aln, optional: true
    file "${sample}/${sample}.phyloFlash.json", emit: summary, optional: true
    path "*.{log,err}"                        , emit: logs, optional: true
    path ".command.*"                         , emit: nf_logs
    path "versions.yml"                       , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def read_opts = meta.single_end ? "-read1 ${reads[0]}" : "-read1 ${reads[0]} -read2 ${reads[1]}"
    """
    mkdir $prefix
    phyloFlash.pl \\
        $options.args \\
        $read_opts \\
        -lib $prefix \\
        -dbhome . \\
        -CPUs $task.cpus

    jsonify-phyloflash.py ${prefix}.phyloFlash > ${prefix}.phyloFlash.json
    mv ${prefix}.* $prefix


    if phyloflash-summary.py ${prefix}/ | grep -q -c "WARNING: Multiple SSUs were assembled by SPAdes"; then
        MULTI="1"
    fi

    if [ "${params.allow_multiple_16s}" == "true" ]; then
        MULTI="0"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloFlash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
