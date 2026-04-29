/**
 * Predict serotypes of Neisseria meningitidis from genome assemblies.
 *
 * This subworkflow uses [meningotype](https://github.com/MDU-PHL/meningotype) to perform
 * in silico serotyping, finetyping and Bexsero antigen sequence typing of *Neisseria meningitidis*
 * strains from assembled genomes. It processes each sample individually and aggregates the
 * results into a single consolidated report.
 *
 * @status stable
 * @keywords neisseria meningitidis, serotype, finetype, bexsero, meningococcal
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation meningotype
 *
 * @modules csvtk_concat, meningotype
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited meningotype results with serogroup, PorA, and FetA predictions
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT                      } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                       } from 'plugin/nf-bactopia'

workflow MENINGOTYPE {
    take:
    assembly: Channel<Record>

    main:
    ch_meningotype = MENINGOTYPE_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_meningotype, 'tsv', [name: 'meningotype']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_meningotype
    run_outputs = ch_csvtk_concat
}
