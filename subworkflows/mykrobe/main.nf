/**
 * Predict antibiotic resistance from sequence reads.
 *
 * This subworkflow uses [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) to predict antibiotic
 * resistance directly from sequencing reads. It provides rapid genotype-based resistance predictions
 * for specific bacterial species.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords bacteria, reads, antimicrobial resistance, genotype prediction
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent, aggregation
 * @citation mykrobe
 *
 * @modules mykrobe_predict, csvtk_concat
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input mykrobe_species
 * Target bacterial species for resistance prediction (e.g., "staphylococcus_aureus",
 * "mycobacterium_tuberculosis", or "enterococcus_faecium").
 *
 * @output csv        Detailed resistance predictions in CSV format for each sample
 * @output json       Machine-readable resistance predictions in JSON format for each sample
 * @output merged_csv Combined resistance predictions from all samples in a single CSV file
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MYKROBE_PREDICT }  from '../../modules/mykrobe/predict/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow MYKROBE {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    mykrobe_species: String

    main:
    MYKROBE_PREDICT(reads, mykrobe_species)
    CSVTK_CONCAT(gather(MYKROBE_PREDICT.out.csv, 'mykrobe'), 'csv', 'csv')

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Set<Path>>> = MYKROBE_PREDICT.out.csv
    json: Channel<Tuple<Map, Set<Path>>> = MYKROBE_PREDICT.out.json
    merged_csv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.csv,
        MYKROBE_PREDICT.out.json,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
