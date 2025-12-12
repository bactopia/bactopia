/**
 * Rapid Haemophilus parasuis serotyping.
 *
 * This subworkflow performs serotyping of Haemophilus parasuis using
 * [HpsuisSero](https://github.com/jimmyliu1326/HpsuisSero), which identifies
 * serotype-specific markers in genome assemblies. The tool provides rapid
 * classification of H. parasuis isolates into their respective serotypes,
 * which is important for epidemiological surveillance and vaccine development.
 *
 * @status stable
 * @keywords Haemophilus, parasuis, serotype, epidemiology
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation hpsuissero
 *
 * @modules csvtk_concat, hpsuissero
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. parasuis serotype prediction
 *
 * @output tsv         HpsuisSero serotype prediction results in TSV format
 * @output merged_tsv  Combined TSV file containing serotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow HPSUISSERO {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    HPSUISSERO_MODULE(assembly)
    CSVTK_CONCAT(gather(HPSUISSERO_MODULE.out.tsv, 'hpsuissero'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = HPSUISSERO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        HPSUISSERO_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
