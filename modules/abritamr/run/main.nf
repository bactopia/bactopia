process ABRITAMR_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.summary_matches.txt")  , emit: matches
    tuple val(meta), path("${prefix}.summary_partials.txt") , emit: partials
    tuple val(meta), path("${prefix}.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("${prefix}.amrfinder.out")        , emit: amrfinder
    tuple val(meta), path("${prefix}.abritamr.txt")         , emit: summary, optional: true
    tuple val(meta), path("*.{log,err}")                    , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                 , emit: nf_begin
    tuple val(meta), path(".command.err")                   , emit: nf_err
    tuple val(meta), path(".command.log")                   , emit: nf_log
    tuple val(meta), path(".command.out")                   , emit: nf_out
    tuple val(meta), path(".command.run")                   , emit: nf_run
    tuple val(meta), path(".command.sh")                    , emit: nf_sh
    tuple val(meta), path(".command.trace")                 , emit: nf_trace
    tuple val(meta), path("versions.yml")                   , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    abritamr run \\
        --contigs $fasta_name \\
        --prefix $prefix \\
        $task.ext.args \\
        --jobs $task.cpus

    # Rename output files to prevent name collisions
    mv ${prefix}/summary_matches.txt ./${prefix}.summary_matches.txt
    mv ${prefix}/summary_partials.txt ./${prefix}.summary_partials.txt
    mv ${prefix}/summary_virulence.txt ./${prefix}.summary_virulence.txt
    mv ${prefix}/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f results/abritamr.txt ]; then
        # This file is not always present
        mv ${prefix}/abritamr.txt ./${prefix}.abritamr.txt
    fi

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
