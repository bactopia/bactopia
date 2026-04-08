/**
 * Identify SCCmec elements in Staphylococcus aureus genomes.
 *
 * This subworkflow uses [SCCmec](https://github.com/rpetit3/sccmec) to identify the
 * Staphylococcal Cassette Chromosome mec (SCCmec) element in *Staphylococcus aureus*
 * assemblies. It predicts the type based on the presence of specific *mec* and *ccr*
 * gene complexes, generating detailed BLAST results and typing information.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation sccmec
 *
 * @modules sccmec, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Main results file with SCCmec typing
 * - `targets`: BLAST results for target sequences
 * - `target_details`: Detailed results for target matches
 * - `regions`: BLAST results for SCCmec regions
 * - `regions_details`: Detailed results for SCCmec region matches
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gatherCsvtk             } from 'plugin/nf-bactopia'

workflow SCCMEC {
    take:
    assembly: Channel<Record>

    main:
    ch_sccmec = SCCMEC_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_sccmec, 'tsv', [name: 'sccmec']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_sccmec
    run_outputs = ch_csvtk_concat
}
