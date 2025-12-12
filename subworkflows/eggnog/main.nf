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
 * @output hits            Best hits from sequence similarity searches
 * @output seed_orthologs  Seed ortholog assignments for each query
 * @output annotations     Detailed functional annotations with GO terms, pathways, and COG categories
 * @output xlsx            Complete annotation results in Excel format
 * @output orthologs       Ortholog assignments in TSV format
 * @output genepred        Gene prediction output file
 * @output gff             Annotation file in GFF3 format
 * @output no_anno         Sequences with no significant annotation
 * @output pfam            PFAM domain annotations
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../modules/eggnog/mapper/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow EGGNOG {
    take:
    proteins: Channel<Tuple<Map, Set<Path>>>
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
    // Individual outputs
    hits: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.hits
    seed_orthologs: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.seed_orthologs
    annotations: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.annotations
    xlsx: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.xlsx
    orthologs: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.orthologs
    genepred: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.genepred
    gff: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.gff
    no_anno: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.no_anno
    pfam: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.pfam

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        EGGNOG_MAPPER.out.hits,
        EGGNOG_MAPPER.out.seed_orthologs,
        EGGNOG_MAPPER.out.annotations,
        EGGNOG_MAPPER.out.xlsx,
        EGGNOG_MAPPER.out.orthologs,
        EGGNOG_MAPPER.out.genepred,
        EGGNOG_MAPPER.out.gff,
        EGGNOG_MAPPER.out.no_anno,
        EGGNOG_MAPPER.out.pfam
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.versions])
}
