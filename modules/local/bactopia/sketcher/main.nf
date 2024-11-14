// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options        = initOptions(params.options ? params.options : [:], 'sketcher')
options.ignore = [".fastq.gz"]
options.btype  = "main"
conda_tools    = "bioconda::bactopia-sketcher=1.0.1"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SKETCHER {
    tag "${meta.id}"
    label "process_low"

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-sketcher:1.0.1--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-sketcher:1.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path mash_db
    path sourmash_db

    output:
    tuple val(meta), path("${prefix}.sig")                        , emit: sig
    tuple val(meta), path("${prefix}-k*.msh")                     , emit: msh
    tuple val(meta), path("${prefix}-mash-refseq88-k21.txt")      , emit: mash
    tuple val(meta), path("${prefix}-sourmash-gtdb-rs207-k31.txt"), emit: sourmash
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = mash_db.getName().endsWith(".xz") ? true : false
    def mash_name = mash_db.getName().replace(".xz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        xz -c -d $mash_db > $mash_name
    fi

    gzip -cd ${fasta} | mash sketch -o ${prefix}-k21 -k 21 ${options.args} -I ${prefix} -
    gzip -cd ${fasta} | mash sketch -o ${prefix}-k31 -k 31 ${options.args} -I ${prefix} -
    sourmash sketch dna ${options.args2} --merge ${prefix} -o ${prefix}.sig ${fasta}

    # Mash Screen
    echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/\t/g' > ${prefix}-mash-refseq88-k21.txt
    gzip -cd ${fasta} | mash screen ${options.args3} -p ${task.cpus} ${mash_name} - | sort -gr >> ${prefix}-mash-refseq88-k21.txt

    # Sourmash classify
    sourmash lca classify --query ${prefix}.sig --db ${sourmash_db} > ${prefix}-sourmash-gtdb-rs207-k31.txt

    # Cleanup
    rm -rf ${mash_name}

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/.*sourmash //;')
    END_VERSIONS
    """
}
