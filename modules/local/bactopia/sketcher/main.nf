// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES      = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options        = initOptions(params.options ? params.options : [:], 'sketcher')
options.ignore = [".fastq.gz"]
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-main/${params.wf}" : params.outdir
conda_tools    = "bioconda::bactopia-sketcher=1.0.0"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SKETCHER {
    tag "${meta.id}"
    label "base_mem_8gb"
    label "minmer_sketch"

    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-sketcher:1.0.0--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-sketcher:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fq)

    output:
    tuple val(meta), path(fq), path("${prefix}.sig"), emit: sketch
    path("${prefix}*.{msh,sig}")
    path("${prefix}.ctx"), optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    fastq = meta.single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    mccortex_fq = meta.single_end ? "-1 ${fq[0]}" : "-2 ${fq[0]}:${fq[1]}"
    m = task.memory.toString().split(' ')[0].toInteger() * 1000 - 500
    """
    gzip -cd ${fastq} | mash sketch -o ${prefix}-k21 -k 21 -s ${params.sketch_size} -r -I ${prefix} -
    gzip -cd ${fastq} | mash sketch -o ${prefix}-k31 -k 31 -s ${params.sketch_size} -r -I ${prefix} -
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=${params.sourmash_scale} --merge ${prefix} -o ${prefix}.sig ${fastq}

    if [[ "${params.count_31mers}" == "true" ]]; then
        mccortex31 build -f -k 31 -s ${prefix} ${mccortex_fq} -t ${task.cpus} -m ${m}mb -q temp_counts
        if [ "${params.keep_singletons}" == "false" ]; then
            # Clean up Cortex file (mostly remove singletons)
            mccortex31 clean -q -B 2 -U2 -T2 -m ${m}mb -o ${prefix}.ctx temp_counts
            rm temp_counts
        else
            mv temp_counts ${prefix}.ctx
        fi
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/sourmash //;')
    END_VERSIONS
    """
}
