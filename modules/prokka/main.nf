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
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input proteins
 * FASTA file of trusted proteins to first annotate from (Optional)
 *
 * @input prodigal_tf
 * Training file to use for gene prediction (Optional)
 *
 * @output record(meta, gff, gbff, fna, faa, ffn, sqn, fsa, tbl, txt, tsv, blastdb, results, logs, nf_logs, versions)
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbff`: Annotation in GenBank format, containing both sequences and annotations
 * - `fna`: Nucleotide FASTA file of the input contig sequences
 * - `faa`: Protein FASTA file of the translated CDS sequences
 * - `ffn`: Nucleotide FASTA file of all prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
 * - `sqn`: An ASN1 format "Sequin" file for submission to GenBank
 * - `fsa`: Nucleotide FASTA file of the input contig sequences, used by tbl2asn
 * - `tbl`: Feature Table file for NCBI submission
 * - `txt`: Summary statistics relating to the annotated features found
 * - `tsv`: Tab-separated file of all features (locus_tag, ftype, len_bp, gene, EC_number, COG, product)
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases of the contigs, genes, and proteins
 */
nextflow.preview.types = true

process PROKKA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Path): Record
    proteins   : Path?
    prodigal_tf: Path?

    stage:
    stageAs "input/*", fna

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        gff: file("${prefix}.{gff,gff.gz}"),
        gbff: file("${prefix}.{gbk,gbk.gz}"),
        fna: file("${prefix}.{fna,fna.gz}"),
        faa: file("${prefix}.{faa,faa.gz}"),
        ffn: file("${prefix}.{ffn,ffn.gz}"),
        sqn: file("${prefix}.{sqn,sqn.gz}"),
        fsa: file("${prefix}.{fsa,fsa.gz}"),
        tbl: file("${prefix}.{tbl,tbl.gz}"),
        txt: file("${prefix}.txt"),
        tsv: file("${prefix}.tsv"),
        blastdb: file("${prefix}-blastdb.tar.gz"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.{gff,gff.gz}"),
            files("${prefix}.{gbk,gbk.gz}"),
            files("${prefix}.{fna,fna.gz}"),
            files("${prefix}.{faa,faa.gz}"),
            files("${prefix}.{ffn,ffn.gz}"),
            files("${prefix}.{sqn,sqn.gz}"),
            files("${prefix}.{fsa,fsa.gz}"),
            files("${prefix}.{tbl,tbl.gz}"),
            files("${prefix}.txt"),
            files("${prefix}.tsv"),
            files("${prefix}-blastdb.tar.gz")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def proteins_opt = proteins != null && proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""
    def prodigal_opt = prodigal_tf != null && prodigal_tf.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigaltf ${prodigal_tf.getName()}" : ""
    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
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
        gzip -c -d ${fna} > ${fna_name}
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
            ${fna_name}
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
            ${fna_name}
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

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi
    rm -rf ${prefix}/ blastdb/ *.pdb *.pjs *.pot *.ptf *.pto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$( echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
