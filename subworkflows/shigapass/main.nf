/**
 * Predict serotypes of Shigella from assemblies.
 *
 * This subworkflow uses [ShigaPass](https://github.com/imanyass/ShigaPass) to predict
 * serotypes of *Shigella* strains from assembled genomes. It analyzes the presence
 * and composition of antigen-encoding genes to classify isolates into their known serotypes.
 *
 * @status stable
 * @keywords shigella, serotype, typing, prediction, antigen genes
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation shigapass
 *
 * @modules shigapass, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: ShigaPass summary results in TSV format
 * - `flex_tsv`: ShigaPass Flex summary results in TSV format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../modules/shigapass/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow SHIGAPASS {
    take:
    assembly: Channel<Record>

    main:
    ch_shigapass = SHIGAPASS_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_shigapass, 'tsv', [name: 'shigapass']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_shigapass
    run_outputs = ch_csvtk_concat
}
