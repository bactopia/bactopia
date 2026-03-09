/**
 * Profile microbial composition using Sylph.
 *
 * This subworkflow estimates microbial composition directly from sequencing reads using
 * [Sylph](https://github.com/xiaoli-dong/sylph). It provides rapid and accurate abundance
 * estimates by comparing k-mer signatures against a reference genome database. Sylph can
 * process both short and long reads, offering taxonomic profiling from species to strain level
 * with confidence estimates for each identification.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenome, profiling, composition, abundance, kmer, taxonomic
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation sylph
 *
 * @modules sylph_profile
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
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

workflow SYLPH {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    database: Path

    main:
    SYLPH_PROFILE(reads, database)

    emit:
    sample_outputs = SYLPH_PROFILE.out
}
