/**
 * Predict phylogroups of Escherichia coli from genome assemblies.
 *
 * This subworkflow uses [ClermontTyping](https://github.com/happykhan/ClermonTyping) to determine
 * the phylogenetic groups of *Escherichia coli* strains from assembled genomes. It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords escherichia coli, phylogroup, typing, clermont
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation clermontyping
 *
 * @modules csvtk_concat, clermontyping
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited ClermonTyping phylogroup assignment results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../modules/clermontyping/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { gather                                } from 'plugin/nf-bactopia'

workflow CLERMONTYPING {
    take:
    assembly: Channel<Record>

    main:
    CLERMONTYPING_MODULE(assembly)
    CSVTK_CONCAT(gather(CLERMONTYPING_MODULE.out, 'clermontyping', field: 'tsv'), 'tsv', 'tsv')

    emit:
    // Per-sample records
    sample_outputs = CLERMONTYPING_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
