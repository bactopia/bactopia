/**
 * In silico prediction of Escherichia coli serotype.
 *
 * This subworkflow performs serotype prediction for Escherichia coli genomes
 * using [ECTyper](https://github.com/phac-nml/ecoli_serotyping), which predicts
 * O and H antigens from whole genome assemblies. The tool identifies specific
 * serotype markers and provides comprehensive serotype classification.
 *
 * @status stable
 * @keywords Escherichia, coli, serotype, O-antigen, H-antigen
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation ectyper
 *
 * @modules csvtk_concat, ectyper
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembly files in FASTA format for E. coli serotype prediction
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited ECTyper serotype prediction results
 * - `txt`: ECTyper detailed results in text format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gatherCsvtk               } from 'plugin/nf-bactopia'

workflow ECTYPER {
    take:
    assembly: Channel<Record>

    main:
    ch_ectyper = ECTYPER_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_ectyper, 'tsv', [name: 'ectyper']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_ectyper
    run_outputs = ch_csvtk_concat
}
