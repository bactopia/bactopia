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
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>
    proteins       : Path?
    prodigal_tf    : Path?

    stage:
    stageAs "input/*", assembly

    output:
    annotations = tuple(meta, files("${prefix}.{fna,fna.gz}"), files("${prefix}.{faa,faa.gz}"), files("${prefix}.{gff,gff.gz}"))
    gff         = tuple(meta, files("${prefix}.{gff,gff.gz}"))
    gbk         = tuple(meta, files("${prefix}.{gbk,gbk.gz}"))
    fna         = tuple(meta, files("${prefix}.{fna,fna.gz}"))
    faa         = tuple(meta, files("${prefix}.{faa,faa.gz}"))
    ffn         = tuple(meta, files("${prefix}.{ffn,ffn.gz}"))
    sqn         = tuple(meta, files("${prefix}.{sqn,sqn.gz}"))
    fsa         = tuple(meta, files("${prefix}.{fsa,fsa.gz}"))
    tbl         = tuple(meta, files("${prefix}.{tbl,tbl.gz}"))
    txt         = tuple(meta, files("${prefix}.txt"))
    tsv         = tuple(meta, files("${prefix}.tsv"))
    blastdb     = tuple(meta, files("${prefix}-blastdb.tar.gz"))
    logs        = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs     = tuple(meta, files(".command.*"))
    versions    = tuple(meta, files("versions.yml"))

    script:
    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""            // TODO: Remove when Path? is fixed
    def prodigal_opt = prodigal_tf?.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigaltf ${prodigal_tf.getName()}" : "" // TODO: Remove when Path? is fixed
    def is_compressed = assembly.toList()[0].getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.toList()[0].getName().replace(".gz", "")
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
