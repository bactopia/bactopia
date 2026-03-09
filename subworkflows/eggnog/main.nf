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
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic, resource-download
 * @citation eggnog_mapper
 *
 * @modules eggnog_download, eggnog_mapper
 *
 * @input tuple(meta, proteins)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `proteins`: Protein sequences in FASTA format for functional annotation
 *
 * @input database
 * Path to pre-downloaded eggNOG database (optional)
 *
 * @input download_eggnog
 * Boolean flag to trigger database download if not provided
 *
 * @output sample_outputs  Per-sample record outputs from EGGNOG_MAPPER
 */
nextflow.preview.types = true

include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../modules/eggnog/mapper/main'

workflow EGGNOG {
    take:
    proteins: Channel<Tuple<Map, Path>>
    database: Path
    download_eggnog: Boolean

    main:
    if (download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(proteins, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(proteins, database)
    }

    emit:
    sample_outputs = EGGNOG_MAPPER.out
}
