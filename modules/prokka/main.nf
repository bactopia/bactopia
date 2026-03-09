/**
 * Annotate prokaryotic genomes.
 *
 * Uses [Prokka](https://github.com/tseemann/prokka) to rapidly annotate bacterial, archaeal,
 * and viral genomes, producing standards-compliant output files including GFF3, GenBank, and Sequin.
 *
 * @status stable
 * @keywords prokka, annotation, prokaryotic, bacteria, genbank, gff
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output,compression,conditional-logic
 * @citation prokka
 *
 * @note Uses EMPTY_* placeholder files for optional parameters
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input proteins
 * FASTA file of trusted proteins to first annotate from (Optional)
 *
 * @input prodigal_tf
 * Training file to use for gene prediction (Optional)
 *
 * @output annotations  A tuple containing the FASTA, protein FASTA, and GFF3 files
 * @output gff          Annotation in GFF3 format, containing both sequences and annotations
 * @output gbk          Annotation in GenBank format, containing both sequences and annotations
 * @output fna          Nucleotide FASTA file of the input contig sequences
 * @output faa          Protein FASTA file of the annotated genes
 * @output ffn          Nucleotide FASTA file of the annotated genes
 * @output sqn          An ASN1 format "Sequin" file for submission to GenBank
 * @output fsa          Nucleotide FASTA file of the input contig sequences (with adjusted headers)
 * @output tbl          Feature table file
 * @output txt          Summary statistics about the annotation
 * @output tsv          Tab-separated file of all features
 * @output blastdb      A compressed tar.gz archive of BLAST databases created from the input
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process PROKKA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record
    proteins    : Path?
    prodigal_tf : Path?

    stage:
    stageAs "input/*", assembly

    output:
    record(
        meta: meta,
        gff: file("${prefix}.{gff,gff.gz}"),
        gbk: file("${prefix}.{gbk,gbk.gz}"),
        fna: file("${prefix}.{fna,fna.gz}"),
        faa: file("${prefix}.{faa,faa.gz}"),
        ffn: file("${prefix}.{ffn,ffn.gz}"),
        sqn: file("${prefix}.{sqn,sqn.gz}"),
        fsa: file("${prefix}.{fsa,fsa.gz}"),
        tbl: file("${prefix}.{tbl,tbl.gz}"),
        txt: file("${prefix}.txt"),
        tsv: file("${prefix}.tsv"),
        blastdb: file("${prefix}-blastdb.tar.gz"),
        results: [
            file("${prefix}.{gff,gff.gz}"),
            file("${prefix}.{gbk,gbk.gz}"),
            file("${prefix}.{fna,fna.gz}"),
            file("${prefix}.{faa,faa.gz}"),
            file("${prefix}.{ffn,ffn.gz}"),
            file("${prefix}.{sqn,sqn.gz}"),
            file("${prefix}.{fsa,fsa.gz}"),
            file("${prefix}.{tbl,tbl.gz}"),
            file("${prefix}.txt"),
            file("${prefix}.tsv"),
            file("${prefix}-blastdb.tar.gz")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""            // TODO: Remove when Path? is fixed
    def prodigal_opt = prodigal_tf?.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigaltf ${prodigal_tf.getName()}" : "" // TODO: Remove when Path? is fixed
    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
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
        gzip -c -d ${assembly} > ${assembly_name}
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
            ${assembly_name}
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
            ${assembly_name}
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
    rm -rf ${assembly_name} ${prefix}/ blastdb/ *.pdb *.pjs *.pot *.ptf *.pto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$( echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
