/**
 * In silico Legionella pneumophila Sequence Based Typing.
 *
 * This subworkflow performs sequence-based typing of Legionella pneumophila
 * using [legsta](https://github.com/tseemann/legsta), which identifies the
 * Sequence Type (ST) based on the seven-locus scheme. The tool analyzes
 * allele profiles and provides epidemiological typing data for outbreak
 * investigation and population studies.
 *
 * @status stable
 * @keywords Legionella, pneumophila, sequence typing, ST, epidemiology
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation legsta
 *
 * @modules csvtk_concat, legsta
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for L. pneumophila sequence typing
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited Legionella pneumophila SBT results with allele numbers and sequence type
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { LEGSTA as LEGSTA_MODULE } from '../../modules/legsta/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gather                  } from 'plugin/nf-bactopia'

workflow LEGSTA {
    take:
    assembly: Channel<Record>

    main:
    LEGSTA_MODULE(assembly)
    CSVTK_CONCAT(gather(LEGSTA_MODULE.out, 'tsv', [name: 'legsta']), 'tsv', 'tsv')

    emit:
    // Per-sample records (contains meta, tsv, results, logs, nf_logs, versions)
    sample_outputs = LEGSTA_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
