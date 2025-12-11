/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat, hicap as hicap_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @input database_dir
 * Channel containing database_dir data
 *
 * @input model_fp
 * Channel containing model_fp data
 *
 * @output gbk        Gbk
 * @output svg        Svg
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow HICAP {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    database_dir: Path?
    model_fp: Path?

    main:
    HICAP_MODULE(fasta, database_dir, model_fp)
    CSVTK_CONCAT(gather(HICAP_MODULE.out.tsv, 'hicap'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gbk: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.gbk
    svg: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.svg
    tsv: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.gbk,
        HICAP_MODULE.out.svg,
        HICAP_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
