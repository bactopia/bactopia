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
 * @output classified     Taxonomically classified reads in FASTA format
 * @output kraken2_report Kraken2 classification report with read counts per taxon
 * @output unclassified   Unclassified reads not assigned to any taxa
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'

workflow KRAKEN2 {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    database: Path

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    sample_outputs = KRAKEN2_MODULE.out
}
