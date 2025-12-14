/**
 * Predict serotypes of Neisseria meningitidis from genome assemblies.
 *
 * This subworkflow uses [meningotype](https://github.com/MDU-PHL/meningotype) to perform
 * in silico serotyping, finetyping and Bexsero antigen sequence typing of *Neisseria meningitidis*
 * strains from assembled genomes. It processes each sample individually and aggregates the
 * results into a single consolidated report.
 *
 * @status stable
 * @keywords neisseria meningitidis, serotype, finetype, bexsero, meningococcal
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation meningotype
 *
 * @modules csvtk_concat, meningotype
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing meningococcal typing results
 * @output merged_tsv  Consolidated TSV file containing meningococcal typing from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT                      } from '../../modules/csvtk/concat/main'
include { flattenPaths                      } from 'plugin/nf-bactopia'
include { gather                            } from 'plugin/nf-bactopia'

workflow MENINGOTYPE {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    MENINGOTYPE_MODULE(assembly)
    CSVTK_CONCAT(gather(MENINGOTYPE_MODULE.out.tsv, 'meningotype'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = MENINGOTYPE_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        MENINGOTYPE_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
