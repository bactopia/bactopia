/**
 * Assign genotypes to Salmonella Typhi genomes.
 *
 * This subworkflow assigns genotypes to Salmonella Typhi genomes based on Mykrobe
 * results using [GenoTyphi](https://github.com/katholt/genotyphi). The workflow first
 * runs Mykrobe for antimicrobial resistance prediction on sequencing reads, then
 * processes the results with GenoTyphi to assign specific genotypes based on
 * the presence of known genetic markers.
 *
 * @status stable
 * @keywords Salmonella, Typhi, genotype, antimicrobial resistance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation genotyphi, mykrobe
 *
 * @modules mykrobe_predict, genotyphi_parse, csvtk_concat
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end sequencing reads for Salmonella Typhi genotype assignment
 *
 * @output tsv         GenoTyphi genotype assignment results in TSV format
 * @output csv         Mykrobe antimicrobial resistance prediction results in CSV format
 * @output json        Mykrobe detailed results in JSON format
 * @output merged_tsv  Combined TSV file containing genotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
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
    tsv: Channel<Tuple<Map, Set<Path>>> = GENOTYPHI_PARSE.out.tsv
    csv: Channel<Tuple<Map, Set<Path>>> = MYKROBE_PREDICT.out.csv
    json: Channel<Tuple<Map, Set<Path>>> = MYKROBE_PREDICT.out.json
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

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
