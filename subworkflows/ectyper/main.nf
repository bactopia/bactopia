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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for E. coli serotype prediction
 *
 * @output sample_outputs  Per-sample records containing meta, tsv, txt, results, logs, nf_logs, versions
 * @output run_outputs   Cross-sample aggregation record
 */
nextflow.preview.types = true

include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gather                    } from 'plugin/nf-bactopia'

workflow ECTYPER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    ECTYPER_MODULE(assembly)
    CSVTK_CONCAT(gather(ECTYPER_MODULE.out, 'ectyper', field: 'tsv'), 'tsv', 'tsv')

    emit:
    // Per-sample records (contains meta, tsv, txt, results, logs, nf_logs, versions)
    sample_outputs = ECTYPER_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
