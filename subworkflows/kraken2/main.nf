/**
 * Classify metagenomic reads using Kraken2.
 *
 * This subworkflow performs taxonomic classification of metagenomic reads using [Kraken2](https://github.com/DerrickWood/kraken2),
 * a fast taxonomic classification system. It assigns taxonomic labels to sequencing reads based on k-mer matching against a reference database.
 *
 * @status stable
 * @keywords metagenomics, taxonomic classification, kraken2, k-mer
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent
 * @citation kraken2
 *
 * @modules kraken2
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Metagenomic reads for taxonomic classification
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
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow KRAKEN2 {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    // Individual outputs
    classified: Channel<Tuple<Map, Set<Path>>> = KRAKEN2_MODULE.out.classified
    kraken2_report: Channel<Tuple<Map, Set<Path>>> = KRAKEN2_MODULE.out.kraken2_report
    unclassified: Channel<Tuple<Map, Set<Path>>> = KRAKEN2_MODULE.out.unclassified

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2_MODULE.out.classified,
        KRAKEN2_MODULE.out.kraken2_report,
        KRAKEN2_MODULE.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2_MODULE.out.versions])
}
