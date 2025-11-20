nextflow.preview.types = true

process PROKKA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    proteins       : Path
    prodigal_tf    : Path

    stage:
    stageAs "input/*", fasta

    output:
    annotations = tuple(meta, file("${prefix}.{fna,fna.gz}"), file("${prefix}.{faa,faa.gz}"), file("${prefix}.{gff,gff.gz}"))
    gff         = tuple(meta, file("${prefix}.{gff,gff.gz}"))
    gbk         = tuple(meta, file("${prefix}.{gbk,gbk.gz}"))
    fna         = tuple(meta, file("${prefix}.{fna,fna.gz}"))
    faa         = tuple(meta, file("${prefix}.{faa,faa.gz}"))
    ffn         = tuple(meta, file("${prefix}.{ffn,ffn.gz}"))
    sqn         = tuple(meta, file("${prefix}.{sqn,sqn.gz}"))
    fsa         = tuple(meta, file("${prefix}.{fsa,fsa.gz}"))
    tbl         = tuple(meta, file("${prefix}.{tbl,tbl.gz}"))
    txt         = tuple(meta, file("${prefix}.txt"))
    tsv         = tuple(meta, file("${prefix}.tsv"))
    blastdb     = tuple(meta, file("${prefix}-blastdb.tar.gz"))
    logs        = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin    = tuple(meta, file(".command.begin"))
    nf_err      = tuple(meta, file(".command.err"))
    nf_log      = tuple(meta, file(".command.log"))
    nf_out      = tuple(meta, file(".command.out"))
    nf_run      = tuple(meta, file(".command.run"))
    nf_sh       = tuple(meta, file(".command.sh"))
    nf_trace    = tuple(meta, file(".command.trace"))
    versions    = tuple(meta, file("versions.yml"))

    script:
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    if (task.ext.wf == "pangenome") {
        meta.scope = "run"
        meta.output_dir = "prokka/${prefix}"
        meta.logs_dir = "prokka/${prefix}/logs"
    }
    else {
        meta.output_dir = "${prefix}/main/annotator/prokka/"
        meta.logs_dir = "${prefix}/main/annotator/prokka/logs/"
    }
    meta.process_name = task.ext.process_name

    // Contig ID must <= 37 characters
    def compliant = task.ext.compliant ? "--compliant" : ""
    def locustag = "--locustag ${prefix}"
    if ("gnl|${task.ext.centre}|${prefix}_100000".length() > 37) {
        locustag = ""
        compliant = "--compliant"
    }
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ "${task.ext.prokka_debug}" == "true" ]; then
        export PROKKA_DBDIR=\$(echo "\$(which prokka | sed "s=/prokka==")/../db")
        env
        mkdir tmp_prokka/
        TMPDIR=tmp_prokka/ bactopia-prokka \\
            ${task.ext.args} \\
            --cpus ${task.cpus} \\
            --prefix ${prefix} \\
            ${compliant} \\
            ${locustag} \\
            ${proteins_opt} \\
            ${prodigal_opt} \\
            ${fasta_name}
        rm -rf tmp_prokka/
    else
        prokka \\
            ${task.ext.args} \\
            --cpus ${task.cpus} \\
            --prefix ${prefix} \\
            ${compliant} \\
            ${locustag} \\
            ${proteins_opt} \\
            ${prodigal_opt} \\
            ${fasta_name}
    fi

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat ${prefix}/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat ${prefix}/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat ${prefix}/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -cvf - blastdb/ | gzip -c > ${prefix}/${prefix}-blastdb.tar.gz

    if [[ "${task.ext.skip_compression}" == "false" ]]; then
        gzip ${prefix}/*.gff
        gzip ${prefix}/*.gbk
        gzip ${prefix}/*.fna
        gzip ${prefix}/*.faa
        gzip ${prefix}/*.ffn
        gzip ${prefix}/*.sqn
        gzip ${prefix}/*.fsa
        gzip ${prefix}/*.tbl
    fi

    mv ${prefix}/* ./

    # Cleanup intermediate files
    rm -rf ${fasta_name} ${prefix}/ blastdb/ *.pdb *.pjs *.pot *.ptf *.pto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$( echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
