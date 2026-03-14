/**
 * Predict spa types of Staphylococcus aureus from genome assemblies.
 *
 * This subworkflow uses [spaTyper](https://github.com/HCGB-IGTP/spaTyper) to predict
 * the spa types of *Staphylococcus aureus* strains from assembled genomes based on
 * the polymorphic X region of the protein A gene (spa). It processes each sample
 * individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords staphylococcus aureus, spa typing, protein a, mrsa
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation spatyper
 *
 * @modules csvtk_concat, spatyper
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input repeats
 * Optional custom repeats database for spa typing
 *
 * @input repeat_order
 * Optional custom repeat order file for spa typing
 *
 * @output sample_outputs
 * - `tsv`: spa typing results in TSV format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SPATYPER as SPATYPER_MODULE } from '../../modules/spatyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gather                      } from 'plugin/nf-bactopia'

workflow SPATYPER {
    take:
    assembly: Channel<Record>
    repeats: Path?
    repeat_order: Path?

    main:
    SPATYPER_MODULE(assembly, repeats, repeat_order)
    CSVTK_CONCAT(gather(SPATYPER_MODULE.out, 'spatyper', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = SPATYPER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
