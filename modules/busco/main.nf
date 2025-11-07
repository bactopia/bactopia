process BUSCO {
    tag "$prefix - $lineage"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)

    output:
    tuple val(meta), path("supplemental/*")       , emit: supplemental
    tuple val(meta), path("${prefix}-summary.txt"), emit: tsv
    tuple val(meta), path("*.{log,err}")          , emit: logs, optional: true
    tuple val(meta), path(".command.begin")       , emit: nf_begin
    tuple val(meta), path(".command.err")         , emit: nf_err
    tuple val(meta), path(".command.log")         , emit: nf_log
    tuple val(meta), path(".command.out")         , emit: nf_out
    tuple val(meta), path(".command.run")         , emit: nf_run
    tuple val(meta), path(".command.sh")          , emit: nf_sh
    tuple val(meta), path(".command.trace")       , emit: nf_trace
    tuple val(meta), path("versions.yml")         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    lineage = task.ext.busco_lineage
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # Have to put FASTA in a directory to force batch mode in busco
    mkdir tmp-fasta
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > tmp-fasta/$fasta_name
    fi

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

    busco \\
        --cpu $task.cpus \\
        --in tmp-fasta/ \\
        --out supplemental \\
        --lineage $lineage \\
        --mode genome \\
        --download_base_url=https://busco-data2.ezlab.org/v5/data \\
        $task.ext.args

    # cleanup output directory structure
    find supplemental/ -name "*.log" | xargs -I {} mv {} ./
    find supplemental/ -type d -path "*logs" | xargs -I {} rm -rf {}
    find supplemental/ -type f -name "*.fna" | xargs -I {} gzip {}
    find supplemental/ -type f -name "*.faa" | xargs -I {} gzip {}
    find supplemental/ -type f -path "*hmmer_output*" -name "*.out" | xargs -I {} gzip {}
    mv supplemental/batch_summary.txt supplemental/${prefix}-summary.txt
    mv supplemental/${fasta_name}/* supplemental/
    rm -rf supplemental/${fasta_name} busco_downloads/ tmp*/

    # Busco outputs additional trailing tabs, clean them up
    sed -i 's/\t\t\t\$//' supplemental/${prefix}-summary.txt
    mv supplemental/${prefix}-summary.txt ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
