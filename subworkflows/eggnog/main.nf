/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules eggnog_mapper, eggnog_download
 *
 * @input faa
 * Channel containing faa data
 *
 * @input database
 * Channel containing database data
 *
 * @input download_eggnog
 * Channel containing download_eggnog data
 *
 * @output hits           Hits
 * @output seed_orthologs Seed Orthologs
 * @output annotations    Annotations
 * @output xlsx           Xlsx
 * @output orthologs      Orthologs
 * @output genepred       Genepred
 * @output gff            Gff
 * @output no_anno        No Anno
 * @output pfam           Pfam
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../modules/eggnog/mapper/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow EGGNOG {
    take:
    faa: Channel<Tuple<Map, Set<Path>>>
    database: Path
    download_eggnog: Boolean

    main:
    if (download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(faa, database)
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
