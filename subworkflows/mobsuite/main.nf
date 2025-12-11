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
 * @modules csvtk_concat, mobsuite_recon
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output txt            Txt
 * @output merged_reports Merged Reports
 * @output chromosome     Chromosome
 * @output contig_report  Contig Report
 * @output plasmids       Plasmids
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow MOBSUITE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    MOBSUITE_RECON(fasta)
    CSVTK_CONCAT(gather(MOBSUITE_RECON.out.txt, 'mobsuite', 'summary'), 'tsv', 'tsv')

    emit:
    txt: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.txt
    merged_reports: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    chromosome: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.chromosome
    contig_report: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.contig_report
    plasmids: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.plasmids

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.txt,
        CSVTK_CONCAT.out.csv,
        MOBSUITE_RECON.out.chromosome,
        MOBSUITE_RECON.out.contig_report,
        MOBSUITE_RECON.out.plasmids
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
