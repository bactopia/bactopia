/**
 * Scripts for finding and processing promoter variants upstream of mcr-1.
 *
 * This subworkflow identifies and characterizes promoter variants upstream of
 * the mcr-1 colistin resistance gene using [mcroni](https://github.com/liampshaw/mcroni).
 * The tool searches for mutations in the promoter region that may affect expression
 * levels of mcr-1, which is important for understanding the regulation of
 * plasmid-mediated colistin resistance.
 *
 * @status stable
 * @keywords mcr-1, colistin, resistance, promoter, variant
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation mcroni
 *
 * @modules csvtk_concat, mcroni
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for mcr-1 promoter analysis
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited mcr-1 gene variation results
 * - `fa`: Extracted mcr-1 gene sequence in FASTA format (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { MCRONI as MCRONI_MODULE } from '../../modules/mcroni/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gather                  } from 'plugin/nf-bactopia'

workflow MCRONI {
    take:
    assembly: Channel<Record>

    main:
    MCRONI_MODULE(assembly)
    CSVTK_CONCAT(gather(MCRONI_MODULE.out, 'mcroni', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = MCRONI_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
