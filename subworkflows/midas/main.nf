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
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of species abundance and coverage
 * - `abundances`: Detailed species abundance profile
 *
 * @output run_outputs
 * - `csv`: Merged species abundance results from all samples
 */
nextflow.preview.types = true

include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT  } from '../../modules/csvtk/concat/main'
include { gather        } from 'plugin/nf-bactopia'

workflow MIDAS {
    take:
    reads: Channel<Record>
    database: Path

    main:
    MIDAS_SPECIES(reads, database)
    CSVTK_CONCAT(gather(MIDAS_SPECIES.out, 'midas', field: 'tsv'), 'tsv', 'tsv')
    emit:
    sample_outputs = MIDAS_SPECIES.out
    run_outputs = CSVTK_CONCAT.out
}
