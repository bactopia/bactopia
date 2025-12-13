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
 * @output hits            Raw search hits (Diamond/MMseqs2) against the eggNOG database
 * @output seed_orthologs  List of identified seed orthologs used for annotation transfer
 * @output annotations     Main tab-delimited annotation file (COGs, KEGG, GO, etc.)
 * @output xlsx            Excel format of the annotations file
 * @output orthologs       List of fine-grained orthologs (optional)
 * @output genepred        Predicted gene sequences (optional)
 * @output gff             Annotations in GFF format (optional)
 * @output no_anno         FASTA file of sequences that failed to be annotated (optional)
 * @output pfam            Raw PFAM domain hits (optional)
 * @output logs            Optional software execution logs containing warnings/errors
 * @output nf_logs         Nextflow execution scripts and logs for debugging
 * @output versions        A YAML formatted file with software versions
 */
nextflow.preview.types = true

process EGGNOG_MAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, proteins) : Tuple<Map, Set<Path>>
    db             : Path

    output:
    hits           = tuple(meta, files("*.emapper.hits"))
    seed_orthologs = tuple(meta, files("*.emapper.seed_orthologs"))
    annotations    = tuple(meta, files("*.emapper.annotations"))
    xlsx           = tuple(meta, files("*.emapper.annotations.xlsx", optional: true))
    orthologs      = tuple(meta, files("*.emapper.orthologs", optional: true))
    genepred       = tuple(meta, files("*.emapper.genepred.fasta", optional: true))
    gff            = tuple(meta, files("*.emapper.gff", optional: true))
    no_anno        = tuple(meta, files("*.emapper.no_annotations.fasta", optional: true))
    pfam           = tuple(meta, files("*.emapper.pfam", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, files("versions.yml"))

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
