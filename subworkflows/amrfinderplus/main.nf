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
 * @modules csvtk_concat, amrfinderplus_run
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @input db
 * Channel containing db data
 *
 * @output report          Report
 * @output merged_tsv      Merged Tsv
 * @output mutation_report Mutation Report
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution logs from all processes
 * @output versions        Aggregated version information from all executed tools
 
 
 
 */
nextflow.preview.types = true

include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT      } from '../../modules/csvtk/concat/main'
include { flattenPaths      } from 'plugin/nf-bactopia'
include { gather            } from 'plugin/nf-bactopia'

workflow AMRFINDERPLUS {
    take:
    fasta: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>>
    db: Path

    main:
    AMRFINDERPLUS_RUN(fasta, db)
    CSVTK_CONCAT(gather(AMRFINDERPLUS_RUN.out.report, 'amrfinderplus'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = AMRFINDERPLUS_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    mutation_report: Channel<Tuple<Map, Path>> = AMRFINDERPLUS_RUN.out.mutation_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.report,
        AMRFINDERPLUS_RUN.out.mutation_report,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
