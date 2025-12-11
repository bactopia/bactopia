/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat, gather as gather_module
 *
 * @input reads
 * Channel containing reads data
 *
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output fastq_only Fastq Only
 * @output raw_fastq  Raw Fastq
 * @output error      Error
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT            } from '../../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow GATHER {
    take:
    reads: Channel<Tuple<Map, Set<Path>, Set<Path>, Path>>

    main:
    GATHER_MODULE(reads)

    // Merge meta values for each sample
    CSVTK_CONCAT(gather(GATHER_MODULE.out.tsv, 'meta'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = GATHER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    fastq_only: Channel<Tuple<Map, Set<Path>>> = GATHER_MODULE.out.fastq_only
    raw_fastq: Channel<Tuple<Map, Set<Path>, Set<Path>>> = GATHER_MODULE.out.raw_fastq
    error: Channel<Tuple<Map, Set<Path>>> = GATHER_MODULE.out.error

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.tsv,
        GATHER_MODULE.out.error,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        GATHER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
