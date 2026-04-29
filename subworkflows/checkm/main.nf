/**
 * Assess metagenome bin completeness using CheckM.
 *
 * This subworkflow evaluates the quality and completeness of metagenome-assembled genomes
 * (MAGs) using [CheckM](https://github.com/Ecogenomics/CheckM). It provides a comprehensive
 * assessment based on lineage-specific marker sets, estimating completeness and contamination
 * of genome bins. The workflow generates detailed reports including marker gene statistics,
 * taxonomy predictions, and quality metrics for each bin.
 *
 * @status stable
 * @keywords metagenome, bin, completeness, contamination, mag, quality
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation checkm
 *
 * @modules csvtk_concat, checkm_lineagewf
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Metagenome-assembled genome bins to evaluate. Each record contains metadata
 *   about the sample and a set of genome bins in FASTA format.
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited CheckM quality assessment results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { CHECKM_LINEAGEWF } from '../../modules/checkm/lineagewf/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gatherCsvtk      } from 'plugin/nf-bactopia'

workflow CHECKM {
    take:
    assembly: Channel<Record>

    main:
    ch_checkm_lineagewf = CHECKM_LINEAGEWF(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_checkm_lineagewf, 'tsv', [name: 'checkm']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_checkm_lineagewf
    run_outputs = ch_csvtk_concat
}
