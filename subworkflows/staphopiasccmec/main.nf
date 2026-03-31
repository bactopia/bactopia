/**
 * Identify SCCmec elements in Staphylococcus aureus genomes using Staphopia method.
 *
 * This subworkflow uses [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec) to
 * identify Staphylococcal Cassette Chromosome mec (SCCmec) elements in *Staphylococcus aureus*
 * assemblies. This is the standalone version of the SCCmec typing method developed for
 * the Staphopia project, which predicts SCCmec types based on the presence of specific
 * *mec* and *ccr* gene complexes.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation staphopiasccmec
 *
 * @modules staphopiasccmec, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: TSV file with SCCmec typing results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE } from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT                              } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                               } from 'plugin/nf-bactopia'

workflow STAPHOPIASCCMEC {
    take:
    assembly: Channel<Record>

    main:
    STAPHOPIASCCMEC_MODULE(assembly)
    CSVTK_CONCAT(gatherCsvtk(STAPHOPIASCCMEC_MODULE.out, 'tsv', [name: 'staphopiasccmec']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = STAPHOPIASCCMEC_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
