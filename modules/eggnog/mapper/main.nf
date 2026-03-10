/**
 * Functional annotation of proteins using eggNOG orthology data.
 *
 * Uses [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) to assign functional annotations
 * to protein sequences. It uses precomputed orthologous groups (OGs) to infer functions like
 * COG categories, KEGG pathways, GO terms, and CAZymes with high precision.
 *
 * @status stable
 * @keywords functional annotation, orthology, cog, kegg, go, proteins, eggnog
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation eggnog_mapper
 *
 * @note Database Required
 * Requires the eggNOG database (including the diamond database and taxonomic data) to be available.
 *
 * @input tuple(meta, proteins)
 * - `meta`: Groovy Map containing sample information
 * - `proteins`: Protein sequences in FASTA format (amino acids)
 *
 * @input db
 * Directory or compressed tarball containing the eggNOG database
 *
 * @output record(meta, hits, seed_orthologs, annotations, xlsx, orthologs, genepred, gff, no_anno, pfam, results, logs, nf_logs, versions)
 * - `hits`: Raw search hits (Diamond/MMseqs2) against the eggNOG database
 * - `seed_orthologs`: List of identified seed orthologs used for annotation transfer
 * - `annotations`: Main tab-delimited annotation file (COGs, KEGG, GO, etc.)
 * - `xlsx`: Excel format of the annotations file
 * - `orthologs`: List of fine-grained orthologs (optional)
 * - `genepred`: Predicted gene sequences (optional)
 * - `gff`: Annotations in GFF format (optional)
 * - `no_anno`: FASTA file of sequences that failed to be annotated (optional)
 * - `pfam`: Raw PFAM domain hits (optional)
 */
nextflow.preview.types = true

process EGGNOG_MAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, proteins: Path): Record
    db                : Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        hits: file("${prefix}.emapper.hits"),
        seed_orthologs: file("${prefix}.emapper.seed_orthologs"),
        annotations: file("${prefix}.emapper.annotations"),
        xlsx: file("${prefix}.emapper.annotations.xlsx", optional: true),
        orthologs: file("${prefix}.emapper.orthologs", optional: true),
        genepred: file("${prefix}.emapper.genepred.fasta", optional: true),
        gff: file("${prefix}.emapper.gff", optional: true),
        no_anno: file("${prefix}.emapper.no_annotations.fasta", optional: true),
        pfam: file("${prefix}.emapper.pfam", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.emapper.hits"),
            files("${prefix}.emapper.seed_orthologs"),
            files("${prefix}.emapper.annotations"),
            files("${prefix}.emapper.annotations.xlsx", optional: true),
            files("${prefix}.emapper.orthologs", optional: true),
            files("${prefix}.emapper.genepred.fasta", optional: true),
            files("${prefix}.emapper.gff", optional: true),
            files("${prefix}.emapper.no_annotations.fasta", optional: true),
            files("${prefix}.emapper.pfam", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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

    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        EGGNOG_DB=\$(find database/ -name "eggnog.db" | sed 's=eggnog.db==')
    else
        EGGNOG_DB=\$(find ${db}/ -name "eggnog.db" | sed 's=eggnog.db==')
    fi

    emapper.py \\
        ${task.ext.args} \\
        --cpu ${task.cpus} \\
        --data_dir \$EGGNOG_DB \\
        --output ${prefix} \\
        -i ${proteins}

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
