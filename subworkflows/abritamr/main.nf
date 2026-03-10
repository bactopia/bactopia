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
 * @input tuple(meta, assembly)
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
include { gather       } from 'plugin/nf-bactopia'

workflow ABRITAMR {
    take:
    assembly: Channel<Record>

    main:
    ABRITAMR_RUN(assembly)
    CSVTK_CONCAT(gather(ABRITAMR_RUN.out, 'abritamr', field: 'summary'), 'tsv', 'tsv')

    emit:
    // Per-sample records (contains meta, summary, matches, partials, virulence, amrfinder, results, logs, nf_logs, versions)
    sample_outputs = ABRITAMR_RUN.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
