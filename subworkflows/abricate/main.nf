/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of [Abricate](https://github.com/tseemann/abricate)
 * to screen genome assemblies for antimicrobial resistance and virulence genes, followed by
 * aggregating the results into a single summary report.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, virulence, workflow
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules abricate_run, abricate_summary
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `report`: A tab-delimited report of hits, for full details please see [Abricate - Output](https://github.com/tseemann/abricate#output)
 *
 * @output run_outputs
 * - `report`: Aggregated tab-delimited summary of Abricate results from all samples
 */
nextflow.preview.types = true

include { ABRICATE_RUN     } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'
include { gather           } from 'plugin/nf-bactopia'

workflow ABRICATE {
    take:
    assembly: Channel<Record>

    main:
    ABRICATE_RUN(assembly)
    ABRICATE_SUMMARY(gather(ABRICATE_RUN.out, 'report', [name: 'abricate']))

    emit:
    // Published outputs
    sample_outputs = ABRICATE_RUN.out
    run_outputs = ABRICATE_SUMMARY.out
}
