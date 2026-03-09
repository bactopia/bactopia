/**
 * Estimate species abundance from metagenomic reads.
 *
 * This subworkflow performs taxonomic classification and abundance estimation using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken). It processes metagenomic reads, classifies them against a reference database,
 * and generates abundance estimates at different taxonomic levels with optional abundance correction.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomic classification, abundance estimation, kraken2, bracken
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, conditional-logic, aggregation
 * @citation kraken2, bracken
 *
 * @modules bracken, csvtk_concat as csvtk_concat_tsv, csvtk_concat as csvtk_concat_adjusted
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
 * @output sample_outputs  Per-sample record outputs from BRACKEN_MODULE
 * @output run_outputs   Combined abundance results across all samples as records
 */
nextflow.preview.types = true

include { BRACKEN as BRACKEN_MODULE             } from '../../modules/bracken/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV      } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../modules/csvtk/concat/main'
include { gather                                } from 'plugin/nf-bactopia'

workflow BRACKEN {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    database: Path

    main:
    BRACKEN_MODULE(reads, database)

    // Merge Bracken Primary/Secondary Species abundance
    CSVTK_CONCAT_TSV(gather(BRACKEN_MODULE.out, 'bracken-species-abundance', field: 'tsv'), 'tsv', 'tsv')

    // Merge Bracken adjusted abundance
    CSVTK_CONCAT_ADJUSTED(gather(BRACKEN_MODULE.out, 'bracken-adjusted', field: 'adjusted_abundances'), 'tsv', 'tsv')

    emit:
    sample_outputs = BRACKEN_MODULE.out
    run_outputs = CSVTK_CONCAT_TSV.out.mix(CSVTK_CONCAT_ADJUSTED.out)
}
