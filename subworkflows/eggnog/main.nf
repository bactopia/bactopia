/**
 * Functional annotation through orthology assignment.
 *
 * This subworkflow performs genome-wide functional annotation using
 * [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), which provides fast
 * functional annotation through orthology assignment. It can optionally download
 * the eggNOG database if not provided. The pipeline predicts orthologs, functional
 * categories, and various annotation formats including GFF, Excel, and detailed reports.
 *
 * @status stable
 * @keywords functional, annotation, orthology, eggNOG, protein domains
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic,resource-download
 * @citation eggnog_mapper
 *
 * @modules eggnog_download, eggnog_mapper
 *
 * @input record(meta, proteins)
 * - `meta`: Groovy Record containing sample information
 * - `proteins`: Protein sequences in FASTA format for functional annotation
 *
 * @input database
 * Path to pre-downloaded eggNOG database (optional)
 *
 * @input download_eggnog
 * Boolean flag to trigger database download if not provided
 *
 * @input save_as_tarball
 * Boolean flag to save downloaded database as tarball
 *
 * @output sample_outputs
 * - `hits`: Raw search hits (Diamond/MMseqs2) against the eggNOG database
 * - `seed_orthologs`: List of identified seed orthologs used for annotation transfer
 * - `annotations`: Main tab-delimited annotation file (COGs, KEGG, GO, etc.)
 * - `xlsx`: Excel format of the annotations file
 * - `orthologs`: List of fine-grained orthologs (optional)
 * - `genepred`: Predicted gene sequences (optional)
 * - `gff`: Annotations in GFF format (optional)
 * - `no_anno`: FASTA file of sequences that failed to be annotated (optional)
 * - `pfam`: Raw PFAM domain hits (optional)
 *
 * @output run_outputs
 */
nextflow.enable.types = true

include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../modules/eggnog/mapper/main'

workflow EGGNOG {
    take:
    proteins: Channel<Record>
    database: Path?
    download_eggnog: Boolean
    save_as_tarball: Boolean

    main:
    ch_eggnog_mapper = channel.empty()
    if (download_eggnog) {
        ch_eggnog_download = EGGNOG_DOWNLOAD()

        if (save_as_tarball) {
            ch_eggnog_mapper = EGGNOG_MAPPER(proteins, ch_eggnog_download.map { r -> r.db_tarball })
        } else {
            ch_eggnog_mapper = EGGNOG_MAPPER(proteins, ch_eggnog_download.map { r -> r.db })
        }
    } else {
        ch_eggnog_mapper = EGGNOG_MAPPER(proteins, database)
    }

    emit:
    // Published outputs
    sample_outputs = ch_eggnog_mapper
    run_outputs = channel.empty()
}
