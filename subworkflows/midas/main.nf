/**
 * Species-level profiling from metagenomic data.
 *
 * This subworkflow estimates strain-level genomic variation from metagenomic data
 * using [MIDAS](https://github.com/snayfach/MIDAS). The pipeline identifies bacterial
 * species abundances and provides strain-level profiling including SNP analysis.
 * It uses a comprehensive reference database for accurate species identification
 * and quantification in complex microbial communities.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, species, profiling, abundance, strain
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation midas
 *
 * @modules csvtk_concat, midas_species, midas_download
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (not supported by MIDAS)
 *
 * @input database
 * MIDAS reference database for species identification and quantification
 *
 * @input download_midas
 * Boolean flag to automatically download the MIDAS database if not available
 *
 * @input save_as_tarball
 * Boolean flag to save downloaded database as tarball
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of species abundance and coverage
 * - `abundances`: Detailed species abundance profile
 *
 * @output run_outputs
 * - `csv`: Merged species abundance results from all samples
 */
nextflow.preview.types = true

include { MIDAS_DOWNLOAD } from '../../modules/midas/download/main'
include { MIDAS_SPECIES  } from '../../modules/midas/species/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { filterWithData } from 'plugin/nf-bactopia'
include { gatherCsvtk    } from 'plugin/nf-bactopia'

workflow MIDAS {
    take:
    reads: Channel<Record>
    database: Path?
    download_midas: Boolean
    save_as_tarball: Boolean

    main:
    def filtered_reads = filterWithData(reads, ['r1', 'r2', 'se'])

    if (download_midas) {
        MIDAS_DOWNLOAD()

        if (save_as_tarball) {
            MIDAS_SPECIES(filtered_reads, MIDAS_DOWNLOAD.out.map { r -> r.db_tarball })
        } else {
            MIDAS_SPECIES(filtered_reads, MIDAS_DOWNLOAD.out.map { r -> r.db })
        }
    } else {
        MIDAS_SPECIES(filtered_reads, database)
    }

    CSVTK_CONCAT(gatherCsvtk(MIDAS_SPECIES.out, 'tsv', [name: 'midas']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = MIDAS_SPECIES.out
    run_outputs = CSVTK_CONCAT.out
}
