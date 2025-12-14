/**
 * Determine multilocus sequence types (MLST) from bacterial assemblies.
 *
 * This subworkflow uses [mlst](https://github.com/tseemann/mlst) to scan assembled
 * contigs against PubMLST typing schemes and determine sequence types (STs). It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords mlst, sequence typing, pubmlst, bacteria
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation mlst
 *
 * @modules csvtk_concat, mlst
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * PubMLST database to use for MLST typing
 *
 * @output tsv         Per-sample TSV files containing MLST typing results
 * @output merged_tsv  Consolidated TSV file containing MLST typing from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT        } from '../../modules/csvtk/concat/main'
include { flattenPaths        } from 'plugin/nf-bactopia'
include { gather              } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    assembly: Channel<Tuple<Map, Path>>
    db: Path

    main:
    MLST_MODULE(assembly, db)
    CSVTK_CONCAT(gather(MLST_MODULE.out.tsv, 'mlst'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = MLST_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        MLST_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
