/**
 * Calculate Mash distances between sequences and a reference.
 *
 * This subworkflow uses [Mash](https://github.com/marbl/Mash) to calculate MinHash-based
 * distances between query sequences and a reference sequence. It creates Mash sketches
 * of the input sequences and computes distance values, then aggregates all distance
 * calculations into a single consolidated report.
 *
 * @status stable
 * @keywords mash, distance, minhash, comparison, reference
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation mash
 *
 * @modules csvtk_concat, mash_dist
 *
 * @input record(meta, seqs)
 * - `meta`: Groovy Map containing sample information
 * - `seqs`: Sequences in FASTA format to compare against reference
 *
 * @input reference
 * Reference sequence in FASTA format for distance calculations
 *
 * @output sample_outputs
 * - `dist`: A tab-delimited summary of the Mash distances and p-values
 *
 * @output run_outputs
 * - `csv`: Merged Mash distance results from all samples
 */
nextflow.preview.types = true

include { MASH_DIST    } from '../../modules/mash/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gather       } from 'plugin/nf-bactopia'

workflow MASHDIST {
    take:
    seqs: Channel<Record>
    reference: Path

    main:
    MASH_DIST(seqs, reference)
    CSVTK_CONCAT(gather(MASH_DIST.out, 'dist', [name: 'mashdist']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = MASH_DIST.out
    run_outputs = CSVTK_CONCAT.out
}
