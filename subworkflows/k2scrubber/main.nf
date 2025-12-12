/**
 * Remove human reads from metagenomic data using Kraken2.
 *
 * This subworkflow identifies and removes human sequences from metagenomic reads using [Kraken2](https://github.com/DerrickWood/kraken2)
 * with a specialized human genome database. It downloads the k2_HPRC human reference database,
 * classifies reads, and separates human (contaminant) sequences from microbial reads for downstream analysis.
 *
 * @status stable
 * @keywords metagenomics, human decontamination, read filtering, kraken2
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, resource-download
 * @citation kraken2
 *
 * @modules wget, kraken2
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Metagenomic reads potentially contaminated with human sequences
 *
 * @output human                Reads classified as human (contaminant sequences)
 * @output kraken2_report       Kraken2 classification report showing human vs. non-human read counts
 * @output scrub_report         Human contamination screening report with detailed statistics
 * @output scrub_special_report Extended human screening report with additional metrics
 * @output scrubbed             Clean metagenomic reads after human sequence removal
 * @output scrubbed_extra       Additional cleaned reads from extended filtering
 * @output results              Aggregated results channel containing all output files
 * @output logs                 Aggregated logs channel containing all execution logs
 * @output nf_logs              Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions             Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { WGET         } from '../../modules/wget/main'
include { KRAKEN2      } from '../../modules/kraken2/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow K2SCRUBBER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    WGET([
        "name": "k2scrubber",
        "save_as": "k2_HPRC_20230810.tar.gz",
        "url": "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1"
    ])
    KRAKEN2(reads, WGET.out.download)

    emit:
    // Individual outputs
    human: Channel<Tuple<Map, Path>> = KRAKEN2.out.classified
    kraken2_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.kraken2_report
    scrub_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_special_report
    scrubbed: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified
    scrubbed_extra: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified_extra

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2.out.classified,
        KRAKEN2.out.kraken2_report,
        KRAKEN2.out.scrub_report,
        KRAKEN2.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([KRAKEN2.out.versions])
}
