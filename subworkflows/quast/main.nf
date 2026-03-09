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
 * @input tuple(meta, assembly, meta_file)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format (Path)
 * - `meta_file`: Meta file containing reference size information (Path)
 *
 * @output sample_outputs  Per-sample record outputs from QUAST_MODULE
 * @output run_outputs   Combined QUAST results across all samples as a record
 */
nextflow.preview.types = true

include { QUAST as QUAST_MODULE } from '../../modules/quast/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow QUAST {
    take:
    fasta: Channel<Tuple<Map, Path, Path>>

    main:
    QUAST_MODULE(fasta)
    CSVTK_CONCAT(gather(QUAST_MODULE.out, 'quast', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = QUAST_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
