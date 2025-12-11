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
 * @modules plasmidfinder as plasmidfinder_module, csvtk_concat
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output tsv         Tsv
 * @output merged_tsv  Merged Tsv
 * @output json        Json
 * @output txt         Txt
 * @output genome_seq  Genome Seq
 * @output plasmid_seq Plasmid Seq
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution logs from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../modules/plasmidfinder/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { flattenPaths                          } from 'plugin/nf-bactopia'
include { gather                                } from 'plugin/nf-bactopia'

workflow PLASMIDFINDER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    PLASMIDFINDER_MODULE(fasta)
    CSVTK_CONCAT(gather(PLASMIDFINDER_MODULE.out.tsv, 'plasmidfinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = PLASMIDFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    json: Channel<Tuple<Map, Path>> = PLASMIDFINDER_MODULE.out.json
    txt: Channel<Tuple<Map, Path>> = PLASMIDFINDER_MODULE.out.txt
    genome_seq: Channel<Tuple<Map, Path>> = PLASMIDFINDER_MODULE.out.genome_seq
    plasmid_seq: Channel<Tuple<Map, Path>> = PLASMIDFINDER_MODULE.out.plasmid_seq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PLASMIDFINDER_MODULE.out.json,
        PLASMIDFINDER_MODULE.out.txt,
        PLASMIDFINDER_MODULE.out.genome_seq,
        PLASMIDFINDER_MODULE.out.plasmid_seq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
