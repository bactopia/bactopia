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
 * @modules mykrobe_predict, genotyphi_parse, csvtk_concat
 *
 * @input reads
 * Channel containing reads data
 *
 * @output tsv        Tsv
 * @output csv        Csv
 * @output json       Json
 * @output merged_tsv Merged Tsv
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { GENOTYPHI_PARSE } from '../../modules/genotyphi/parse/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow GENOTYPHI {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    MYKROBE_PREDICT(reads, "typhi")
    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    CSVTK_CONCAT(gather(GENOTYPHI_PARSE.out.tsv, 'genotyphi'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = GENOTYPHI_PARSE.out.tsv
    csv: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.csv
    json: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.json
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GENOTYPHI_PARSE.out.tsv,
        MYKROBE_PREDICT.out.csv,
        MYKROBE_PREDICT.out.json,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.logs,
        GENOTYPHI_PARSE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GENOTYPHI_PARSE.out.nf_logs,
        MYKROBE_PREDICT.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.versions,
        GENOTYPHI_PARSE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
