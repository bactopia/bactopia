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
 * @modules sistr as sistr_module, csvtk_concat
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output tsv          Tsv
 * @output merged_tsv   Merged Tsv
 * @output allele_fasta Allele Fasta
 * @output allele_json  Allele Json
 * @output cgmlst_csv   Cgmlst Csv
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution logs from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SISTR as SISTR_MODULE } from '../../modules/sistr/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow SISTR {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    SISTR_MODULE(fasta)
    CSVTK_CONCAT(gather(SISTR_MODULE.out.tsv, 'sistr'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SISTR_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    allele_fasta: Channel<Tuple<Map, Path>> = SISTR_MODULE.out.allele_fasta
    allele_json: Channel<Tuple<Map, Path>> = SISTR_MODULE.out.allele_json
    cgmlst_csv: Channel<Tuple<Map, Path>> = SISTR_MODULE.out.cgmlst_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        SISTR_MODULE.out.allele_fasta,
        SISTR_MODULE.out.allele_json,
        SISTR_MODULE.out.cgmlst_csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
