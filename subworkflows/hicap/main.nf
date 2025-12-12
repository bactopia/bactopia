/**
 * In silico serotyping of the Haemophilus influenzae capsule locus.
 *
 * This subworkflow performs serotyping of Haemophilus influenzae by analyzing
 * the capsule locus using [hicap](https://github.com/scwatts/hicap). The tool
 * identifies capsule type and produces detailed annotations of the cap locus,
 * including graphical representation in SVG format and annotation in GenBank format.
 *
 * @status stable
 * @keywords Haemophilus, influenzae, serotype, capsule, cap locus
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-input, aggregation
 * @citation hicap
 *
 * @modules csvtk_concat, hicap
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. influenzae serotype prediction
 *
 * @input database_dir
 * Directory containing hicap reference database files (optional)
 *
 * @input model_fp
 * Path to HMM model file for improved detection (optional)
 *
 * @output gbk         Annotated capsule locus in GenBank format
 * @output svg         Graphical representation of the capsule locus structure
 * @output tsv         hicap serotype prediction results in TSV format
 * @output merged_tsv  Combined TSV file containing serotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow HICAP {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    database_dir: Path?
    model_fp: Path?

    main:
    HICAP_MODULE(assembly, database_dir, model_fp)
    CSVTK_CONCAT(gather(HICAP_MODULE.out.tsv, 'hicap'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gbk: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.gbk
    svg: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.svg
    tsv: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.gbk,
        HICAP_MODULE.out.svg,
        HICAP_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
