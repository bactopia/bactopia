/**
 * Calculate Average Nucleotide Identity (ANI) between genomes.
 *
 * This subworkflow uses [FastANI](https://github.com/ParBLiSS/FastANI) to compute
 * whole-genome Average Nucleotide Identity (ANI) values between query genomes
 * and reference genomes. ANI is a robust measure of genomic similarity used for
 * species delineation in microbial taxonomy. The results are aggregated into
 * a single consolidated report.
 *
 * @status stable
 * @keywords ani, average nucleotide identity, taxonomy, species, comparison
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation fastani
 *
 * @modules csvtk_concat, fastani
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Query genomes in FASTA format for ANI calculation
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Reference genomes in FASTA format for ANI calculation
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `tsv`: A tab-delimited summary of the ANI scores, matched fragments, and total fragments
 * - `csv`: Aggregated results in CSV format
 */
// bactopia-lint: ignore S015
nextflow.preview.types = true

include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gatherFields              } from 'plugin/nf-bactopia'
include { gatherCsvtk               } from 'plugin/nf-bactopia'

workflow FASTANI {
    take:
    query: Channel<Record>
    reference: Channel<Record>

    main:
    ch_ref = reference.map { r -> r.fna }
    FASTANI_MODULE(gatherFields(query, [fna: 'query'], [name: 'fastani']), ch_ref)
    CSVTK_CONCAT(gatherCsvtk(FASTANI_MODULE.out, 'tsv', [name: 'fastani']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = FASTANI_MODULE.out.mix(CSVTK_CONCAT.out)
}
