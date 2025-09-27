process CLONALFRAMEML {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(msa), path(newick)

    output:
    tuple val(meta), path("${task.ext.process_name}/*.emsim.txt")                      , emit: emsim, optional: true
    tuple val(meta), path("${task.ext.process_name}/*.em.txt")                         , emit: em
    tuple val(meta), path("${task.ext.process_name}/*.importation_status.txt")         , emit: status
    tuple val(meta), path("${task.ext.process_name}/*.labelled_tree.newick")           , emit: newick
    tuple val(meta), path("${task.ext.process_name}/*.ML_sequence.fasta.gz")           , emit: fasta
    tuple val(meta), path("${task.ext.process_name}/*.position_cross_reference.txt.gz"), emit: pos_ref
    tuple val(meta), path("*.masked.aln.gz"), emit: masked_aln
    tuple val(meta), path("*.{log,err}")    , emit: logs, optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process_name}"
    meta.name = prefix
    meta.output_dir = "${task.ext.rundir}/"
    meta.logs_dir = "${task.ext.rundir}/${task.ext.process_name}/logs/"
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

    # Organize 
    mkdir ${task.ext.process_name}/
    mv *.emsim.txt ${task.ext.process_name}/
    mv *.em.txt ${task.ext.process_name}/
    mv *.importation_status.txt ${task.ext.process_name}/
    mv *.labelled_tree.newick ${task.ext.process_name}/
    mv *.ML_sequence.fasta.gz ${task.ext.process_name}/
    mv *.position_cross_reference.txt.gz ${task.ext.process_name}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clonalframeml: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
        maskrc-svg: \$( echo \$(maskrc-svg.py --version 2>&1) | sed 's/^.*maskrc-svg.py //' )
    END_VERSIONS
    """
}
