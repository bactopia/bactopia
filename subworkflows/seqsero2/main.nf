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
 * @modules csvtk_concat, seqsero2 as seqsero2_module
 *
 * @input seqs
 * Channel containing seqs data
 *
 * @output tsv        Tsv
 * @output txt        Txt
 * @output merged_tsv Merged Tsv
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SEQSERO2 as SEQSERO2_MODULE } from '../../modules/seqsero2/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SEQSERO2 {
    take:
    seqs: Channel<Tuple<Map, Set<Path>>>

    main:
    SEQSERO2_MODULE(seqs)
    CSVTK_CONCAT(gather(SEQSERO2_MODULE.out.tsv, 'seqsero2'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SEQSERO2_MODULE.out.tsv
    txt: Channel<Tuple<Map, Path>> = SEQSERO2_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.tsv,
        SEQSERO2_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
