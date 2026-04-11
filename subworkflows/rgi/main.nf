/**
 * Predict antimicrobial resistance from protein or nucleotide data.
 *
 * This subworkflow uses the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) to predict
 * resistomes based on homology and SNP models. It includes analysis of resistance genes,
 * creation of summary visualizations, and aggregation of results across samples.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, resistome, homology
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation rgi
 *
 * @modules rgi_main, rgi_heatmap, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format for resistome prediction
 *
 * @output sample_outputs
 * - `tsv`: RGI results in tab-separated format
 * - `json`: RGI results in JSON format (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
// bactopia-lint: ignore S015
nextflow.preview.types = true

include { RGI_MAIN     } from '../../modules/rgi/main/main'
include { RGI_HEATMAP  } from '../../modules/rgi/heatmap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gather       } from 'plugin/nf-bactopia'
include { gatherCsvtk  } from 'plugin/nf-bactopia'

workflow RGI {
    take:
    assembly: Channel<Record>

    main:
    ch_rgi_main = RGI_MAIN(assembly)
    ch_rgi_heatmap = RGI_HEATMAP(gather(ch_rgi_main, 'json', [name: 'rgi']))
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_rgi_main, 'tsv', [name: 'rgi']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_rgi_main
    run_outputs = ch_rgi_heatmap.mix(ch_csvtk_concat)
}
