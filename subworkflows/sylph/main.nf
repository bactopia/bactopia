/**
 * Profile microbial composition using Sylph.
 *
 * This subworkflow estimates microbial composition directly from sequencing reads using
 * [Sylph](https://github.com/xiaoli-dong/sylph). It provides rapid and accurate abundance
 * estimates by comparing k-mer signatures against a reference genome database. Sylph can
 * process both short and long reads, offering taxonomic profiling from species to strain level
 * with confidence estimates for each identification.
 *
 * @status stable
 * @keywords metagenome, profiling, composition, abundance, kmer, taxonomic
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation sylph
 *
 * @modules sylph_profile
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Sequencing reads to profile. Each tuple contains a set of paired-end or single-end reads in FASTQ format
 *
 * @input database
 * Path to Sylph reference database directory containing pre-computed
 * k-mer signatures of reference genomes for taxonomic classification.
 *
 * @output tsv      Per-sample microbial composition profile with relative abundances
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow SYLPH {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    SYLPH_PROFILE(reads, database)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SYLPH_PROFILE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([SYLPH_PROFILE.out.tsv])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SYLPH_PROFILE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SYLPH_PROFILE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SYLPH_PROFILE.out.versions])
}
