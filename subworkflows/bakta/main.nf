/**
 * Rapid bacterial genome annotation.
 *
 * This subworkflow uses [Bakta](https://github.com/oschwengers/bakta) to provide
 * rapid, comprehensive annotation of bacterial genomes. It can download and prepare
 * the Bakta database on-demand or use a pre-existing database. The workflow processes
 * each sample individually, producing multiple output formats including GFF3, GenBank,
 * protein sequences, nucleotide sequences, and a BLAST database.
 *
 * @status stable
 * @keywords bacteria, annotation, genome, functional annotation, taxonomy
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent, conditional-logic, archive-output, resource-download
 * @citation bakta
 *
 * @modules bakta_download, bakta_run
 *
 * @note Database can be automatically downloaded or provided as pre-existing tarball
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input database
 * Optional pre-existing Bakta database path
 *
 * @input download_bakta
 * Boolean flag to trigger automatic database download
 *
 * @input save_as_tarball
 * Boolean flag to save downloaded database as tarball
 *
 * @input proteins
 * Optional trusted protein sequences for homology search
 *
 * @input prodigal_tf
 * Optional Prodigal training file for improved gene prediction
 *
 * @input replicons
 * Optional replicon sequences for plasmid identification
 *
 * @output annotations       Complete annotation package containing GFF3, GBK, FAA, FFN, FNA, and other formats
 * @output tsv               Tab-delimited summary of annotation results
 * @output txt               Text summary of annotation results
 * @output embl              EMBL format annotation file
 * @output faa               Amino acid sequences of predicted proteins
 * @output ffn               Nucleotide sequences of predicted features
 * @output fna               Nucleotide sequences of contigs
 * @output gbff              GenBank format annotation file
 * @output gff               Genome annotation in GFF3 format
 * @output hypotheticals_faa Amino acid sequences of hypothetical proteins
 * @output hypotheticals_tsv Tab-delimited summary of hypothetical proteins
 * @output blastdb           BLAST database created from annotation results
 * @output results           Aggregated results channel containing all output files
 * @output logs              Aggregated logs channel containing all execution logs
 * @output nf_logs           Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions          Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN      } from '../../modules/bakta/run/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow BAKTA {
    take:
    assembly: Channel<Tuple<Map, Path>>
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
            BAKTA_RUN(assembly, BAKTA_DOWNLOAD.out.db_tarball, proteins, prodigal_tf, replicons)
        } else {
            BAKTA_RUN(assembly, BAKTA_DOWNLOAD.out.db, proteins, prodigal_tf, replicons)
        }
    } else {
        BAKTA_RUN(assembly, database, proteins, prodigal_tf, replicons)
    }

    emit:
    // Individual outputs
    annotations: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> = BAKTA_RUN.out.annotations
    tsv: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.tsv
    txt: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.txt
    embl: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.embl
    faa: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.faa
    ffn: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.ffn
    fna: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.fna
    gbff: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.gbff
    gff: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.gff
    hypotheticals_faa: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.hypotheticals_tsv
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
