/**
 * Assign genotypes to Salmonella Typhi genomes.
 *
 * This subworkflow assigns genotypes to Salmonella Typhi genomes based on Mykrobe
 * results using [GenoTyphi](https://github.com/katholt/genotyphi). The workflow first
 * runs Mykrobe for antimicrobial resistance prediction on sequencing reads, then
 * processes the results with GenoTyphi to assign specific genotypes based on
 * the presence of known genetic markers.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords Salmonella, Typhi, genotype, antimicrobial resistance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation genotyphi, mykrobe
 *
 * @modules mykrobe_predict, genotyphi_parse, csvtk_concat
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited report containing the assigned GenoTyphi genotype
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
// bactopia-lint: ignore S015
nextflow.preview.types = true

include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { GENOTYPHI_PARSE } from '../../modules/genotyphi/parse/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk     } from 'plugin/nf-bactopia'

workflow GENOTYPHI {
    take:
    reads: Channel<Record>

    main:
    ch_mykrobe_predict = MYKROBE_PREDICT(reads, "typhi")
    ch_genotyphi_parse = GENOTYPHI_PARSE(ch_mykrobe_predict.map { r ->
        record(meta: r.meta, json: r.json)
    })
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_genotyphi_parse, 'tsv', [name: 'genotyphi']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_mykrobe_predict.mix(ch_genotyphi_parse)
    run_outputs = ch_csvtk_concat
}
