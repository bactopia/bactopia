/**
 * Identify antimicrobial resistance genes using AMRFinderPlus.
 *
 * This subworkflow uses [abriTAMR](https://github.com/MDU-PHL/abritamr) to identify
 * antimicrobial resistance genes in bacterial genomes. It runs AMRFinderPlus on
 * each sample and collates the results into functional classes, producing detailed
 * reports on resistance genes, partial matches, and virulence factors.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, amr, amrfinderplus, classification
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation abritamr
 *
 * @modules abritamr_run, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `summary`: Tab-delimited NATA-accredited AMR report summary
 * - `matches`: Tab-delimited list of matched AMR genes
 * - `partials`: Tab-delimited list of partially matched AMR genes
 * - `virulence`: Tab-delimited list of detected virulence genes
 * - `amrfinder`: Raw AMRFinderPlus output
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gatherCsvtk       } from 'plugin/nf-bactopia'

workflow ABRITAMR {
    take:
    assembly: Channel<Record>

    main:
    ABRITAMR_RUN(assembly)
    CSVTK_CONCAT(gatherCsvtk(ABRITAMR_RUN.out, 'summary', [name: 'abritamr']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ABRITAMR_RUN.out
    run_outputs = CSVTK_CONCAT.out
}
