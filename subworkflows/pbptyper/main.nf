/**
 * Predict penicillin binding protein (PBP) types of Streptococcus pneumoniae from genome assemblies.
 *
 * This subworkflow uses [pbptyper](https://github.com/rpetit3/pbptyper) to predict
 * the penicillin binding protein (PBP) types and predict antimicrobial susceptibility
 * of *Streptococcus pneumoniae* strains from assembled genomes. It processes each sample
 * individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, pbp typing, penicillin, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation pbptyper
 *
 * @modules pbptyper, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing PBP typing results
 * @output merged_tsv  Consolidated TSV file containing PBP typing from all samples
 * @output blast       Per-sample BLAST results for PBP sequence matches
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow PBPTYPER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    PBPTYPER_MODULE(assembly)
    CSVTK_CONCAT(gather(PBPTYPER_MODULE.out.tsv, 'pbptyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = PBPTYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    blast: Channel<Tuple<Map, Path>> = PBPTYPER_MODULE.out.blast

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PBPTYPER_MODULE.out.blast
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        PBPTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
