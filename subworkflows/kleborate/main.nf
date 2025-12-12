/**
 * Genotyping tool for Klebsiella pneumoniae and its related species complex.
 *
 * This subworkflow performs comprehensive genotyping of Klebsiella pneumoniae
 * and related species using [Kleborate](https://github.com/katholt/Kleborate). The tool
 * identifies capsular (K) and O-antigen (L) loci, virulence factors, and acquired
 * antimicrobial resistance genes, providing a detailed genotype for surveillance
 * and epidemiological studies.
 *
 * @status stable
 * @keywords Klebsiella, pneumoniae, genotyping, virulence, capsule
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation kleborate
 *
 * @modules csvtk_concat, kleborate
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for Klebsiella genotyping
 *
 * @output tsv         Kleborate genotyping results with K and L loci, virulence, and resistance genes
 * @output merged_tsv  Combined TSV file containing genotyping results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { KLEBORATE as KLEBORATE_MODULE } from '../../modules/kleborate/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow KLEBORATE {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    KLEBORATE_MODULE(assembly)
    CSVTK_CONCAT(gather(KLEBORATE_MODULE.out.txt, 'kleborate', 'txt'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = KLEBORATE_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        KLEBORATE_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
