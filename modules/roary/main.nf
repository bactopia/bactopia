process ROARY {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("roary/*")                        , emit: results
    tuple val(meta), path("core-genome.aln.gz")             , emit: aln, optional: true
    tuple val(meta), path("roary/gene_presence_absence.csv"), emit: csv, optional: true
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
    meta.output_dir = "${task.ext.rundir}/"
    meta.logs_dir = "${task.ext.rundir}/${task.ext.process_name}/logs"
    meta.process_name = task.ext.process_name
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs -r gunzip

    # Roary only supports .gff extension, will need to adjust for gff3 (Bakta) files
    # https://github.com/sanger-pathogens/Roary/blob/master/lib/Bio/Roary/PrepareInputFiles.pm#L82
    # note for later: "xargs -r" will not run if no files are found
    find gff/ -name "*.gff3.gz" | xargs -r gunzip
    find gff/ -name "*.gff3" -print0 | while read -d \$'\0' file; do mv "\$file" "\${file%.gff3}.gff"; done

    roary \\
        ${task.ext.args} \\
        -p $task.cpus \\
        -f supplemental/ \\
        gff/*.gff

    gzip supplemental/*.aln
    gzip supplemental/*.fa

    if [[ -f "supplemental/core_gene_alignment.aln.gz" ]]; then
        cp supplemental/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    # clean up
    rm -rf gff/

    mv supplemental/ roary/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """
}
