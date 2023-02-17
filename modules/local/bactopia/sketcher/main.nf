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
    path mash_db
    path sourmash_db

    output:
    tuple val(meta), path(fq), path("${prefix}.sig")              , emit: sig
    tuple val(meta), path(fq), path("*.msh")                      , emit: msh
    tuple val(meta), path("${prefix}-mash-refseq88-k21.txt")      , emit: mash
    tuple val(meta), path("${prefix}-sourmash-gtdb-rs207-k31.txt"), emit: sourmash
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

    # Mash Screen
    echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/\t/g' > ${prefix}-mash-refseq88-k21.txt
    mash screen ${options.args3} -p ${task.cpus} ${mash_db} ${prefix}-k21.msh | sort -gr >> ${prefix}-mash-refseq88-k21.txt

    # Sourmash classify
    sourmash lca classify --query ${prefix}.sig --db ${sourmash_db} > ${prefix}-sourmash-gtdb-rs207-k31.txt

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/sourmash //;')
    END_VERSIONS
    """
}
