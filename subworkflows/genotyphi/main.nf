/**
 * Assign genotypes to Salmonella Typhi genomes.
 *
 * This subworkflow assigns genotypes to Salmonella Typhi genomes based on Mykrobe
 * results using [GenoTyphi](https://github.com/katholt/genotyphi). The workflow first
 * runs Mykrobe for antimicrobial resistance prediction on sequencing reads, then
 * processes the results with GenoTyphi to assign specific genotypes based on
 * the presence of known genetic markers.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords Salmonella, Typhi, genotype, antimicrobial resistance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation genotyphi, mykrobe
 *
 * @modules mykrobe_predict, genotyphi_parse, csvtk_concat
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited report containing the assigned GenoTyphi genotype
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { GENOTYPHI_PARSE } from '../../modules/genotyphi/parse/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { gather          } from 'plugin/nf-bactopia'

workflow GENOTYPHI {
    take:
    reads: Channel<Record>

    main:
    MYKROBE_PREDICT(reads, "typhi")
    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    CSVTK_CONCAT(gather(GENOTYPHI_PARSE.out, 'tsv', [name: 'genotyphi']), 'tsv', 'tsv')

    // Bridge MYKROBE_PREDICT tuple outputs into records (not yet converted)
    ch_mykrobe_samples = MYKROBE_PREDICT.out.csv
        .join(MYKROBE_PREDICT.out.json)
        .join(MYKROBE_PREDICT.out.logs)
        .join(MYKROBE_PREDICT.out.nf_logs)
        .join(MYKROBE_PREDICT.out.versions)
        .map { meta, csv, json, logs, nf_logs, versions ->
            record(meta: meta, csv: csv, json: json, results: [csv, json], logs: logs, nf_logs: nf_logs, versions: versions)
        }

    emit:
    sample_outputs = GENOTYPHI_PARSE.out.mix(ch_mykrobe_samples)
    run_outputs = CSVTK_CONCAT.out
}
