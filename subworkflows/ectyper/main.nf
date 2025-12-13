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
 * @output tsv         ECTyper serotype prediction results in TSV format
 * @output txt         Additional serotype prediction output in text format
 * @output merged_tsv  Combined TSV file containing serotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow ECTYPER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    ECTYPER_MODULE(assembly)
    CSVTK_CONCAT(gather(ECTYPER_MODULE.out.tsv, 'ectyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = ECTYPER_MODULE.out.tsv
    txt: Channel<Tuple<Map, Set<Path>>> = ECTYPER_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.tsv,
        ECTYPER_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        ECTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
