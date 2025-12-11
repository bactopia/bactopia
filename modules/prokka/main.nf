/**
 * Annotate prokaryotic genomes.
 *
 * This process executes prokka to perform analysis
 *
 * @status stable
 * @keywords prokka, annotation, prokaryotic
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, path-workarounds
 * @citation prokka
 *
 * @note Uses EMPTY_* placeholder files for optional parameters
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Input FASTA file (contigs or complete genome)
 *
 * @input proteins
 * Fasta file of trusted proteins to first annotate from
 *
 * @input prodigal_tf
 * Training file to use for gene prediction
 *
 * @output annotations Annotations
 * @output gff         annotation in GFF3 format, containing both sequences and annotations
 * @output gbk         annotation in GenBank format, containing both sequences and annotations
 * @output fna         nucleotide FASTA file of the input contig sequences
 * @output faa         protein FASTA file of the annotated genes
 * @output ffn         nucleotide FASTA file of the annotated genes
 * @output sqn         an ASN1 format "Sequin" file for submission to GenBank
 * @output fsa         nucleotide FASTA file of the input contig sequences (with adjusted sequence headers)
 * @output tbl         feature table file
 * @output txt         summary statistics about the annotation
 * @output tsv         tab-separated file of all features (locus_tag,ftype,gene,EC_number,COG,product)
 * @output blastdb     a compressed tar.gz archive of BLAST databases created from the input sequences and annotations
 * @output logs        Optional tool execution logs
 * @output nf_logs     Nextflow execution logs
 * @output versions    Software version information (YAML format)
 */
nextflow.preview.types = true

process PROKKA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>
    proteins       : Path?
    prodigal_tf    : Path?

    stage:
    stageAs "input/*", fasta

    output:
    annotations = tuple(meta, files("${prefix}.{fna,fna.gz}"), files("${prefix}.{faa,faa.gz}"), files("${prefix}.{gff,gff.gz}"))
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
    blastdb     = tuple(meta, files("${prefix}-blastdb.tar.gz"))
    logs        = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs     = tuple(meta, files(".command.*"))
    versions    = tuple(meta, file("versions.yml"))

    script:
    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""            // TODO: Remove when Path? is fixed
    def prodigal_opt = prodigal_tf.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigaltf ${prodigal_tf.getName()}" : "" // TODO: Remove when Path? is fixed
    def is_compressed = fasta.toList()[0].getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.toList()[0].getName().replace(".gz", "")
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
    echo "${proteins}"
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
