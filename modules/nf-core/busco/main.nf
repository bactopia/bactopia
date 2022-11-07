// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'busco')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::busco=5.3.2"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BUSCO {
    tag "$lineage"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir: lineage) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0' :
        'quay.io/biocontainers/busco:5.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path('tmp_input/*')
    each lineage

    output:
    tuple val(meta), path("${lineage}/")                      , emit: results
    tuple val(meta), path("${lineage}/${lineage}-summary.txt"), emit: tsv
    path "*.{log,err}"                                        , emit: logs, optional: true
    path ".command.*"                                         , emit: nf_logs
    path "versions.yml"                                       , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .fna.gz )
        else
            cp "\$FASTA" \$( basename "\$FASTA" .fna )
        fi
    done
    cd ..

    busco \\
        --cpu $task.cpus \\
        --in "\$INPUT_SEQS" \\
        --out $lineage \\
        --lineage $lineage \\
        --mode genome \\
        --download_base_url=https://busco-data2.ezlab.org/v5/data \\
        $options.args

    mv ${lineage}/batch_summary.txt ${lineage}/${lineage}-summary.txt

    # Busco outputs additional trailing tabs, clean them up
    sed -i 's/\t\t\t\$//' ${lineage}/${lineage}-summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
