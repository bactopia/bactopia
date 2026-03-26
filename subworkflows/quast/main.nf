/**
 * Evaluate assembly quality using QUAST.
 *
 * This subworkflow assesses genome assembly quality using [QUAST](https://quast.sourceforge.net/)
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
 * @input record(meta, assembly, meta_file)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format (Path)
 * - `meta_file`: Meta file containing reference size information (Path)
 *
 * @output sample_outputs
 * - `tsv`: Transposed report in TSV format
 * - `supplemental`: Supplemental files including plots and HTML reports
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { QUAST as QUAST_MODULE } from '../../modules/quast/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow QUAST {
    take:
    fasta: Channel<Record>

    main:
    QUAST_MODULE(fasta)
    CSVTK_CONCAT(gather(QUAST_MODULE.out, 'tsv', [name: 'quast']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = QUAST_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
