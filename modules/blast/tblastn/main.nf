nextflow.preview.types = true

process BLAST_TBLASTN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, blastdb) : Tuple<Map, Path>
    query            : Path

    output:
    tsv      = tuple(meta, file('*.tblastn.tsv'))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
    versions = tuple(meta, file("versions.yml"))

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
    // genes -> ffn, contigs -> fna
    def db_type = task.ext.use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf ${blastdb}
    
    ${which_cat} ${query} | \\
    tblastn \\
        -num_threads ${task.cpus} \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        ${task.ext.args} \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "${outcols}" | sed 's/<TAB>/\t/g' > ${prefix}.tblastn.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastn.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
