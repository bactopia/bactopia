/**
 * Perform multi-antigen sequence typing of Neisseria gonorrhoeae from genome assemblies.
 *
 * This subworkflow uses [ngmaster](https://github.com/MDU-PHL/ngmaster) to perform
 * in silico multi-antigen sequence typing (NG-MAST) for *Neisseria gonorrhoeae*
 * strains from assembled genomes. It processes each sample individually and aggregates
 * the results into a single consolidated report.
 *
 * @status stable
 * @keywords neisseria gonorrhoeae, ng-mast, typing, gonococcal, antigen
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation ngmaster
 *
 * @modules csvtk_concat, ngmaster
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing NG-MAST typing results
 * @output merged_tsv  Consolidated TSV file containing NG-MAST typing from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { NGMASTER as NGMASTER_MODULE } from '../../modules/ngmaster/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow NGMASTER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    NGMASTER_MODULE(assembly)
    CSVTK_CONCAT(gather(NGMASTER_MODULE.out.tsv, 'ngmaster'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = NGMASTER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        NGMASTER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        NGMASTER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        NGMASTER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        NGMASTER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
