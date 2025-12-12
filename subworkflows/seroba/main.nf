/**
 * k-mer based pipeline to identify the serotype of Streptococcus pneumoniae.
 *
 * This subworkflow performs serotyping of Streptococcus pneumoniae from Illumina
 * next-generation sequencing reads using [Seroba](https://github.com/sanger-pathogens/seroba).
 * The tool uses a k-mer based approach to rapidly classify pneumococcal isolates into
 * their respective serotypes based on the capsular polysaccharide synthesis locus.
 *
 * @status stable
 * @keywords Streptococcus, pneumoniae, serotype, k-mer, capsular
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation seroba
 *
 * @modules csvtk_concat, seroba_run
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for S. pneumoniae serotype prediction
 *
 * @output tsv         Seroba serotype prediction results in TSV format
 * @output txt         Detailed serotype assignment report
 * @output merged_tsv  Combined TSV file containing serotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SEROBA_RUN   } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow SEROBA {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    SEROBA_RUN(assembly)
    CSVTK_CONCAT(gather(SEROBA_RUN.out.tsv, 'seroba'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SEROBA_RUN.out.tsv
    txt: Channel<Tuple<Map, Path>> = SEROBA_RUN.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.tsv,
        SEROBA_RUN.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
