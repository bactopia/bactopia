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
 * @input record(meta, r1, r2, se, lr)
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
 * @output sample_outputs
 * - `csv`: AMR predictions in machine-readable CSV format
 * - `json`: Detailed AMR prediction results in JSON format
 *
 * @output run_outputs
 * - `csv`: Merged AMR predictions from all samples
 */
nextflow.preview.types = true

include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { gather         } from 'plugin/nf-bactopia'

workflow MYKROBE {
    take:
    reads: Channel<Record>
    mykrobe_species: String

    main:
    MYKROBE_PREDICT(reads, mykrobe_species)
    CSVTK_CONCAT(gather(MYKROBE_PREDICT.out, 'csv', [name: 'mykrobe']), 'csv', 'csv')
    emit:
    sample_outputs = MYKROBE_PREDICT.out
    run_outputs = CSVTK_CONCAT.out
}
