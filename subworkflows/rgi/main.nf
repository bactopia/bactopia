/**
 * Predict antimicrobial resistance from protein or nucleotide data.
 *
 * This subworkflow uses the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) to predict
 * resistomes based on homology and SNP models. It includes analysis of resistance genes,
 * creation of summary visualizations, and aggregation of results across samples.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, resistome, homology
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation rgi
 *
 * @modules rgi_main, rgi_heatmap, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format for resistome prediction
 *
 * @output sample_outputs
 * - `tsv`: RGI results in tab-separated format
 * - `json`: RGI results in JSON format (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { RGI_MAIN     } from '../../modules/rgi/main/main'
include { RGI_HEATMAP  } from '../../modules/rgi/heatmap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gather       } from 'plugin/nf-bactopia'

workflow RGI {
    take:
    assembly: Channel<Record>

    main:
    RGI_MAIN(assembly)
    CSVTK_CONCAT(gather(RGI_MAIN.out, 'rgi', field: 'tsv'), 'tsv', 'tsv')
    RGI_HEATMAP(gather(RGI_MAIN.out, 'rgi', field: 'json'))

    emit:
    sample_outputs = RGI_MAIN.out
    run_outputs = CSVTK_CONCAT.out
}
