nextflow.preview.types = true

process CLONALFRAMEML {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, msa, newick) : Tuple<Map, Path, Path>

    output:
    emsim      = tuple(meta, file("${task.ext.process_name}/*.emsim.txt", optional: true))
    em         = tuple(meta, file("${task.ext.process_name}/*.em.txt"))
    status     = tuple(meta, file("${task.ext.process_name}/*.importation_status.txt"))
    newick     = tuple(meta, file("${task.ext.process_name}/*.labelled_tree.newick"))
    fasta      = tuple(meta, file("${task.ext.process_name}/*.ML_sequence.fasta.gz"))
    pos_ref    = tuple(meta, file("${task.ext.process_name}/*.position_cross_reference.txt.gz"))
    masked_aln = tuple(meta, file("*.masked.aln.gz"))
    logs       = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin   = tuple(meta, file(".command.begin"))
    nf_err     = tuple(meta, file(".command.err"))
    nf_log     = tuple(meta, file(".command.log"))
    nf_out     = tuple(meta, file(".command.out"))
    nf_run     = tuple(meta, file(".command.run"))
    nf_sh      = tuple(meta, file(".command.sh"))
    nf_trace   = tuple(meta, file(".command.trace"))
    versions   = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.ext.process_name}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${task.ext.process_name}/logs/"
    meta.process_name = task.ext.process_name
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${msa} > ${msa_name}
    fi

    ClonalFrameML \\
        ${newick} \\
        ${msa_name} \\
        ${prefix} \\
        ${task.ext.args}

    maskrc-svg.py ${prefix} --aln ${msa_name} --symbol '-' --out ${prefix}.masked.aln
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
