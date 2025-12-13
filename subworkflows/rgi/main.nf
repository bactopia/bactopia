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
 * @output tsv        RGI prediction results in TSV format for each sample
 * @output merged_tsv  Combined RGI results from all samples in a single TSV file
 * @output json       Detailed RGI predictions in JSON format for each sample
 * @output heatmap    Heatmap visualization of resistance gene presence across all samples
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { RGI_MAIN     } from '../../modules/rgi/main/main'
include { RGI_HEATMAP  } from '../../modules/rgi/heatmap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow RGI {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    RGI_MAIN(assembly)
    CSVTK_CONCAT(gather(RGI_MAIN.out.tsv, 'rgi'), 'tsv', 'tsv')
    RGI_HEATMAP(gather(RGI_MAIN.out.json, 'rgi', 'json'))

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = RGI_MAIN.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    json: Channel<Tuple<Map, Set<Path>>> = RGI_MAIN.out.json
    heatmap: Channel<Tuple<Map, Set<Path>>> = RGI_HEATMAP.out.heatmap

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.tsv,
        RGI_HEATMAP.out.heatmap,
        CSVTK_CONCAT.out.csv,
        RGI_MAIN.out.json
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.logs,
        RGI_HEATMAP.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.nf_logs,
        RGI_HEATMAP.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.versions,
        RGI_HEATMAP.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
