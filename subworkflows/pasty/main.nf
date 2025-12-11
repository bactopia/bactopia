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
 * @modules csvtk_concat, pasty as pasty_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output blast      Blast
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow PASTY {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    PASTY_MODULE(fasta)
    CSVTK_CONCAT(gather(PASTY_MODULE.out.tsv, 'pasty'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = PASTY_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    blast: Channel<Tuple<Map, Path>> = PASTY_MODULE.out.blast

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.tsv,
        PASTY_MODULE.out.blast,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
