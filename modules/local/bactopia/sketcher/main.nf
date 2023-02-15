// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES      = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options        = initOptions(params.options ? params.options : [:], 'sketcher')
options.ignore = [".fastq.gz"]
options.btype  = options.btype ?: "main"
conda_tools    = "bioconda::bactopia-sketcher=1.0.0"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SKETCHER {
    tag "${meta.id}"
    label "base_mem_8gb"
    label "minmer_sketch"

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-sketcher:1.0.0--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-sketcher:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fq)

    output:
    tuple val(meta), path(fq), path("${prefix}.sig"), emit: sig
    path("${prefix}*.{msh,sig}")
    path("${prefix}.ctx"), optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    fastq = meta.single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    """
    gzip -cd ${fastq} | mash sketch -o ${prefix}-k21 -k 21 ${options.args} -r -I ${prefix} -
    gzip -cd ${fastq} | mash sketch -o ${prefix}-k31 -k 31 ${options.args} -r -I ${prefix} -
    sourmash sketch dna ${options.args2} --merge ${prefix} -o ${prefix}.sig ${fastq}

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/sourmash //;')
    END_VERSIONS
    """
}
