process CLONALFRAMEML {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(msa), path(newick)

    output:
    tuple val(meta), path("*.emsim.txt")                      , emit: emsim, optional: true
    tuple val(meta), path("*.em.txt")                         , emit: em
    tuple val(meta), path("*.importation_status.txt")         , emit: status
    tuple val(meta), path("*.labelled_tree.newick")           , emit: newick
    tuple val(meta), path("*.ML_sequence.fasta.gz")           , emit: fasta
    tuple val(meta), path("*.position_cross_reference.txt.gz"), emit: pos_ref
    tuple val(meta), path("*.masked.aln.gz")                  , emit: masked_aln
    tuple val(meta), path("*.{log,err}")                      , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                   , emit: nf_begin
    tuple val(meta), path(".command.err")                     , emit: nf_err
    tuple val(meta), path(".command.log")                     , emit: nf_log
    tuple val(meta), path(".command.out")                     , emit: nf_out
    tuple val(meta), path(".command.run")                     , emit: nf_run
    tuple val(meta), path(".command.sh")                      , emit: nf_sh
    tuple val(meta), path(".command.trace")                   , emit: nf_trace
    tuple val(meta), path("versions.yml")                     , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    ClonalFrameML \\
        $newick \\
        $msa_name \\
        $prefix \\
        $task.ext.args

    maskrc-svg.py $prefix --aln ${msa_name} --symbol '-' --out ${prefix}.masked.aln
    gzip ${prefix}.masked.aln

    # Cleanup
    rm ${msa_name}
    gzip *.ML_sequence.fasta
    gzip *.position_cross_reference.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clonalframeml: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
        maskrc-svg: \$( echo \$(maskrc-svg.py --version 2>&1) | sed 's/^.*maskrc-svg.py //' )
    END_VERSIONS
    """
}
