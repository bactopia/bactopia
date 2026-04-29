/**
 * Predict phylogroups of Escherichia coli from genome assemblies.
 *
 * This subworkflow uses [ClermontTyping](https://github.com/happykhan/ClermonTyping) to determine
 * the phylogenetic groups of *Escherichia coli* strains from assembled genomes. It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords escherichia coli, phylogroup, typing, clermont
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation clermontyping
 *
 * @modules csvtk_concat, clermontyping
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited ClermonTyping phylogroup assignment results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../modules/clermontyping/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                           } from 'plugin/nf-bactopia'

workflow CLERMONTYPING {
    take:
    assembly: Channel<Record>

    main:
    ch_clermontyping = CLERMONTYPING_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_clermontyping, 'tsv', [name: 'clermontyping']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_clermontyping
    run_outputs = ch_csvtk_concat
}
