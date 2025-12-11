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
 * @modules bakta_download, bakta_run
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @input database
 * Channel containing database data
 *
 * @input download_bakta
 * Channel containing download_bakta data
 *
 * @input save_as_tarball
 * Channel containing save_as_tarball data
 *
 * @input proteins
 * Channel containing proteins data
 *
 * @input prodigal_tf
 * Channel containing prodigal_tf data
 *
 * @input replicons
 * Channel containing replicons data
 *
 * @output annotations       Annotations
 * @output tsv               Tsv
 * @output txt               Txt
 * @output embl              Embl
 * @output faa               Faa
 * @output ffn               Ffn
 * @output fna               Fna
 * @output gbff              Gbff
 * @output gff               Gff
 * @output hypotheticals_faa Hypotheticals Faa
 * @output hypotheticals_tsv Hypotheticals Tsv
 * @output blastdb           Blastdb
 * @output results           Aggregated results channel containing all output files
 * @output logs              Aggregated logs channel containing all execution logs
 * @output nf_logs           Aggregated Nextflow execution logs from all processes
 * @output versions          Aggregated version information from all executed tools
 
 */
nextflow.preview.types = true

include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN      } from '../../modules/bakta/run/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow BAKTA {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    database: Path?
    download_bakta: Boolean
    save_as_tarball: Boolean
    proteins: Path?
    prodigal_tf: Path?
    replicons: Path?

    main:
    if (download_bakta) {
        // Force BAKTA_DOWNLOAD to wait
        BAKTA_DOWNLOAD()

        if (save_as_tarball) {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db_tarball, proteins, prodigal_tf, replicons)
        } else {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db, proteins, prodigal_tf, replicons)
        }
    } else {
        BAKTA_RUN(fasta, database, proteins, prodigal_tf, replicons)
    }

    emit:
    // Individual outputs
    annotations: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> = BAKTA_RUN.out.annotations
    tsv: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.tsv
    txt: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.txt
    embl: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.embl
    faa: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.faa
    ffn: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.ffn
    fna: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.fna
    gbff: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.gbff
    gff: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.gff
    hypotheticals_faa: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.hypotheticals_tsv
    blastdb: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.blastdb

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BAKTA_RUN.out.tsv,
        BAKTA_RUN.out.txt,
        BAKTA_RUN.out.embl,
        BAKTA_RUN.out.faa,
        BAKTA_RUN.out.ffn,
        BAKTA_RUN.out.fna,
        BAKTA_RUN.out.gbff,
        BAKTA_RUN.out.gff,
        BAKTA_RUN.out.hypotheticals_faa,
        BAKTA_RUN.out.hypotheticals_tsv,
        BAKTA_RUN.out.blastdb
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.versions])
}
