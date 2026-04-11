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
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent,conditional-logic,archive-output,resource-download
 * @citation bakta
 *
 * @modules bakta_download, bakta_run
 *
 * @note Database can be automatically downloaded or provided as pre-existing tarball
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
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
 * @output sample_outputs
 * - `embl`: Annotations and sequences in EMBL format
 * - `faa`: CDS/sORF amino acid sequences as FASTA
 * - `ffn`: Feature nucleotide sequences as FASTA
 * - `fna`: Replicon/contig DNA sequences as FASTA
 * - `gbff`: Annotations and sequences in GenBank format
 * - `gff`: Annotations and sequences in GFF3 format
 * - `hypotheticals_tsv`: Further information on hypothetical protein CDS as tab-separated values
 * - `hypotheticals_faa`: Hypothetical protein CDS amino acid sequences as FASTA
 * - `tsv`: Annotations as simple human readable tab-separated values
 * - `txt`: Broad summary of Bakta annotations
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases of the contigs, genes, and proteins
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN      } from '../../modules/bakta/run/main'
include { filterWithData } from 'plugin/nf-bactopia'
include { gatherCsvtk    } from 'plugin/nf-bactopia'

workflow BAKTA {
    take:
    assembly: Channel<Record>
    database: Value<Path?>
    download_bakta: Boolean
    save_as_tarball: Boolean
    proteins: Value<Path?>
    prodigal_tf: Value<Path?>
    replicons: Value<Path?>

    main:
    ch_bakta_run = channel.empty()
    if (download_bakta) {
        // Force BAKTA_DOWNLOAD to wait
        ch_bakta_download = BAKTA_DOWNLOAD()

        if (save_as_tarball) {
            ch_bakta_run = BAKTA_RUN(assembly, ch_bakta_download.map { r -> r.db_tarball }, proteins, prodigal_tf, replicons)
        } else {
            ch_bakta_run = BAKTA_RUN(assembly, ch_bakta_download.map { r -> r.db }, proteins, prodigal_tf, replicons)
        }
    } else {
        ch_bakta_run = BAKTA_RUN(assembly, database, proteins, prodigal_tf, replicons)
    }

    emit:
    // Downstream inputs
    annotations = filterWithData(ch_bakta_run, ['fna', 'faa', 'gff'])

    // Published outputs
    sample_outputs = ch_bakta_run
    run_outputs = channel.empty()
}
