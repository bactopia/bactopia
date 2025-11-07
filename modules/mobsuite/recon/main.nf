process MOBSUITE_RECON {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}-chromosome.fasta.gz"), emit: chromosome
    tuple val(meta), path("${prefix}-contig_report.txt")  , emit: contig_report
    tuple val(meta), path("plasmid_*.fasta.gz")           , emit: plasmids, optional: true
    tuple val(meta), path("${prefix}-mobtyper.txt")       , emit: txt, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        ${task.ext.args} \\
        --num_threads $task.cpus \\
        --outdir supplemental \\
        --sample_id $prefix

    if [[ -f "supplemental/mobtyper_results.txt" ]]; then
        mv supplemental/mobtyper_results.txt ${prefix}-mobtyper.txt
    fi

    if [[ -f "supplemental/chromosome.fasta" ]]; then
        mv supplemental/chromosome.fasta ${prefix}-chromosome.fasta
        gzip ${prefix}-chromosome.fasta
    fi

    if [[ -f "supplemental/contig_report.txt" ]]; then
        mv supplemental/contig_report.txt ${prefix}-contig_report.txt
    fi

    # Cleanup
    gzip supplemental/*.fasta
    mv supplemental/*.fasta.gz ./
    rm -rf ${fasta_name} supplemental/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
