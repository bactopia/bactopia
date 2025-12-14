/**
 * Species-level profiling from metagenomic data.
 *
 * This subworkflow estimates strain-level genomic variation from metagenomic data
 * using [MIDAS](https://github.com/snayfach/MIDAS). The pipeline identifies bacterial
 * species abundances and provides strain-level profiling including SNP analysis.
 * It uses a comprehensive reference database for accurate species identification
 * and quantification in complex microbial communities.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, species, profiling, abundance, strain
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation midas
 *
 * @modules csvtk_concat, midas_species
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (not supported by MIDAS)
 *
 * @input database
 * MIDAS reference database for species identification and quantification
 *
 * @output tsv         Species identification and abundance results in TSV format
 * @output merged_tsv  Combined TSV file containing species results from all samples
 * @output abundances  Detailed abundance profiles for detected species
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT  } from '../../modules/csvtk/concat/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow MIDAS {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    database: Path

    main:
    MIDAS_SPECIES(reads, database)
    CSVTK_CONCAT(gather(MIDAS_SPECIES.out.tsv, 'midas'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = MIDAS_SPECIES.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    abundances: Channel<Tuple<Map, Set<Path>>> = MIDAS_SPECIES.out.abundances

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.tsv,
        MIDAS_SPECIES.out.abundances,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
