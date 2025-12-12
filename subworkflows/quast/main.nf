/**
 * Evaluate assembly quality using QUAST.
 *
 * This subworkflow assesses genome assembly quality using [QUAST](http://quast.sourceforge.net/)
 * (Quality Assessment Tool for Genome Assemblies). It provides comprehensive metrics
 * including N50, L50, GC content, total length, and other quality statistics. The workflow
 * generates both individual sample reports and a combined summary for comparative analysis
 * across all assemblies.
 *
 * @status stable
 * @keywords assembly, quality, assessment, metrics, n50, evaluation
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation quast
 *
 * @modules csvtk_concat, quast
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Genome assemblies to evaluate. Each tuple contains metadata about the sample
 *   and sets of assembled contigs in FASTA format (primary and optional secondary assemblies).
 *
 * @output tsv        Per-sample QUAST quality assessment results in TSV format
 * @output merged_tsv Combined QUAST results summary across all samples
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { QUAST as QUAST_MODULE } from '../../modules/quast/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow QUAST {
    take:
    fasta: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    QUAST_MODULE(fasta)
    CSVTK_CONCAT(gather(QUAST_MODULE.out.tsv, 'quast'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = QUAST_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.supplemental,
        QUAST_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
