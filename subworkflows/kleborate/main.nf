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
 * @output sample_outputs  Per-sample records containing meta, tsv, results, logs, nf_logs, and versions
 * @output run_outputs   Merged record containing meta, csv, results, logs, nf_logs, and versions
 */
nextflow.preview.types = true

include { KLEBORATE as KLEBORATE_MODULE } from '../../modules/kleborate/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gather                        } from 'plugin/nf-bactopia'

workflow KLEBORATE {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    KLEBORATE_MODULE(assembly)
    CSVTK_CONCAT(gather(KLEBORATE_MODULE.out, 'kleborate', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = KLEBORATE_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
