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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation meningotype
 *
 * @modules csvtk_concat, meningotype
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited meningotype results with serogroup, PorA, and FetA predictions
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT                      } from '../../modules/csvtk/concat/main'
include { gather                            } from 'plugin/nf-bactopia'

workflow MENINGOTYPE {
    take:
    assembly: Channel<Record>

    main:
    MENINGOTYPE_MODULE(assembly)
    CSVTK_CONCAT(gather(MENINGOTYPE_MODULE.out, 'meningotype', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = MENINGOTYPE_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
