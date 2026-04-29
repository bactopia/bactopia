/**
 * Predict serotypes of Shigella and EIEC from assemblies.
 *
 * This subworkflow uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to predict
 * serotypes of *Shigella* and Enteroinvasive *E. coli* (EIEC) from assembled genomes.
 * It uses a cluster-informed approach to identify specific serotype markers and classify
 * isolates based on their antigenic profiles.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, typing, cluster analysis
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation shigeifinder
 *
 * @modules shigeifinder, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: ShigEiFinder results in TSV format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../modules/shigeifinder/main'
include { CSVTK_CONCAT                        } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                         } from 'plugin/nf-bactopia'

workflow SHIGEIFINDER {
    take:
    assembly: Channel<Record>

    main:
    ch_shigeifinder = SHIGEIFINDER_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_shigeifinder, 'tsv', [name: 'shigeifinder']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_shigeifinder
    run_outputs = ch_csvtk_concat
}
