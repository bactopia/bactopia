process MIDAS_SPECIES {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}.midas.tsv"), emit: tsv
    tuple val(meta), path("*.abundances.txt")   , emit: abundances
    tuple val(meta), path("*.{log,err}")        , emit: logs, optional: true
    tuple val(meta), path(".command.begin")     , emit: nf_begin
    tuple val(meta), path(".command.err")       , emit: nf_err
    tuple val(meta), path(".command.log")       , emit: nf_log
    tuple val(meta), path(".command.out")       , emit: nf_out
    tuple val(meta), path(".command.run")       , emit: nf_run
    tuple val(meta), path(".command.sh")        , emit: nf_sh
    tuple val(meta), path(".command.trace")     , emit: nf_trace
    tuple val(meta), path("versions.yml")       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def read_opts = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        MIDAS_DB=\$(find database/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    else
        MIDAS_DB=\$(find $db/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    fi

    run_midas.py \\
        species \\
        results \\
        $read_opts \\
        $args \\
        -d \${MIDAS_DB} \\
        -t $task.cpus

    mv results/species/species_profile.txt ${prefix}.midas.abundances.txt
    midas-summary.py ${prefix} ${prefix}.midas.abundances.txt

    # Cleanup
    rm -rf results/
    if [ "$is_tarball" == "true" ]; then
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: $VERSION
    END_VERSIONS
    """
}
