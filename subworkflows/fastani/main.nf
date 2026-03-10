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
 * @input tuple(meta, query)
 * - `meta`: Groovy Map containing sample information
 * - `query`: Query genomes in FASTA format for ANI calculation
 *
 * @input tuple(meta, reference)
 * - `meta`: Groovy Map containing sample information
 * - `reference`: Reference genomes in FASTA format for ANI calculation
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of the ANI scores, matched fragments, and total fragments
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gather                    } from 'plugin/nf-bactopia'

workflow FASTANI {
    take:
    query: Channel<Record>
    reference: Channel<Record>

    main:
    FASTANI_MODULE(gather(query, 'query', 'fasta'), gather(reference, 'reference', 'fasta'))
    CSVTK_CONCAT(gather(FASTANI_MODULE.out, 'fastani', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = FASTANI_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
