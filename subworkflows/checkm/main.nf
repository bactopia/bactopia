/**
 * Assess metagenome bin completeness using CheckM.
 *
 * This subworkflow evaluates the quality and completeness of metagenome-assembled genomes
 * (MAGs) using [CheckM](https://github.com/Ecogenomics/CheckM). It provides a comprehensive
 * assessment based on lineage-specific marker sets, estimating completeness and contamination
 * of genome bins. The workflow generates detailed reports including marker gene statistics,
 * taxonomy predictions, and quality metrics for each bin.
 *
 * @status stable
 * @keywords metagenome, bin, completeness, contamination, mag, quality
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation checkm
 *
 * @modules csvtk_concat, checkm_lineagewf
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Metagenome-assembled genome bins to evaluate. Each tuple contains metadata
 *   about the sample and a set of genome bins in FASTA format.
 *
 * @output report        Per-bin CheckM quality assessment results in TSV format
 * @output merged_reports Combined CheckM results summary across all bins
 * @output results       Aggregated results channel containing all output files
 * @output logs          Aggregated logs channel containing all execution logs
 * @output nf_logs       Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions      Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { CHECKM_LINEAGEWF } from '../../modules/checkm/lineagewf/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow CHECKM {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    CHECKM_LINEAGEWF(assembly)
    CSVTK_CONCAT(gather(CHECKM_LINEAGEWF.out.tsv, 'checkm'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = CHECKM_LINEAGEWF.out.tsv
    merged_reports: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.tsv,
        CHECKM_LINEAGEWF.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
