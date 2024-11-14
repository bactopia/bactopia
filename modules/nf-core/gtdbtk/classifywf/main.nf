// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'gtdb')
options.btype = "tools"
conda_tools   = "bioconda::gtdbtk=2.4.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.4.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fna, stageAs: 'fna-tmp/*')
    path db, stageAs: 'gtdb/*'

    output:
    path "results/*"                                        , emit: results
    tuple val(meta), path("results/${prefix}.*.summary.tsv"), emit: tsv
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        export GTDBTK_DATA_PATH="\$(realpath \$(find database/ -path "*metadata*" -name "metadata.txt" | sed 's=/metadata/metadata.txt=='))"
    else
        export GTDBTK_DATA_PATH="\$(readlink $db)"
    fi
    mkdir fna
    cp -L fna-tmp/* fna/
    find fna/ -name "*.fna.gz" | xargs gunzip

    gtdbtk classify_wf \\
        $options.args \\
        --cpus $task.cpus \\
        --pplacer_cpus $task.cpus \\
        --genome_dir ./fna \\
        --out_dir results \\
        --skip_ani_screen \\
        --prefix ${prefix}
    mv results/*.log ./

    # Cleanup
    if [ "$is_tarball" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi
    if [ "${params.gtdb_keep_msa}" == "false" ]; then
        # Delete MSA of submitted and reference genomes.
        rm -rf results/align/*.msa.fasta.gz
    fi
    rm -rf fna/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
