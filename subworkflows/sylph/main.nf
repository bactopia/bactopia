/**
 * Profile microbial composition using Sylph.
 *
 * This subworkflow estimates microbial composition directly from sequencing reads using
 * [Sylph](https://github.com/xiaoli-dong/sylph). It provides rapid and accurate abundance
 * estimates by comparing k-mer signatures against a reference genome database. Sylph can
 * process both short and long reads, offering taxonomic profiling from species to strain level
 * with confidence estimates for each identification.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenome, profiling, composition, abundance, kmer, taxonomic
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation sylph
 *
 * @modules sylph_profile, csvtk_concat
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to Sylph reference database directory containing pre-computed
 * k-mer signatures of reference genomes for taxonomic classification.
 *
 * @output sample_outputs
 * - `tsv`: TSV file with profiling results
 *
 * @output run_outputs
 * - `csv`: Aggregated profiling results in CSV format
 */
nextflow.preview.types = true

include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'
include { CSVTK_CONCAT  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk   } from 'plugin/nf-bactopia'

workflow SYLPH {
    take:
    reads: Channel<Record>
    database: Value<Path>

    main:
    SYLPH_PROFILE(reads, database)
    CSVTK_CONCAT(gatherCsvtk(SYLPH_PROFILE.out, 'tsv', [name: 'sylph']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = SYLPH_PROFILE.out
    run_outputs = CSVTK_CONCAT.out
}
