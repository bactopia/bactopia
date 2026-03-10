/**
 * Classify metagenomic reads using Kraken2.
 *
 * This subworkflow performs taxonomic classification of metagenomic reads using [Kraken2](https://github.com/DerrickWood/kraken2),
 * a fast taxonomic classification system. It assigns taxonomic labels to sequencing reads based on k-mer matching against a reference database.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomic classification, kraken2, k-mer
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent
 * @citation kraken2
 *
 * @modules kraken2
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to the Kraken2 database for taxonomic classification.
 *
 * @output sample_outputs
 * - `kraken2_report`: Standard Kraken2 report containing taxonomic abundance counts
 * - `scrub_report`: Summary report of reads removed during host scrubbing (optional)
 * - `special_meta`: A simplified metadata map for internal use
 * - `classified`: Reads assigned to a taxon in the database (FASTQ)
 * - `unclassified`: Reads NOT assigned to any taxon (FASTQ)
 * - `classified_extra`: Duplicate classified channel with placeholder for pipeline routing
 * - `unclassified_extra`: Duplicate unclassified channel with placeholder for pipeline routing
 */
nextflow.preview.types = true

include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'

workflow KRAKEN2 {
    take:
    reads: Channel<Record>
    database: Path

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    sample_outputs = KRAKEN2_MODULE.out
}
