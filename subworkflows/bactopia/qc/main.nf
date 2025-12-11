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
 * @modules qc as qc_module
 *
 * @input reads
 * Channel containing reads data
 *
 * @input adapters
 * Channel containing adapters data
 *
 * @input phix
 * Channel containing phix data
 *
 * @output fastq       Fastq
 * @output fastq_only  Fastq Only
 * @output txt         Txt
 * @output error       Error
 * @output error_fastq Error Fastq
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution logs from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { QC as QC_MODULE } from '../../../modules/bactopia/qc/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow QC {
    take:
    reads: Channel<Tuple<Map, Set<Path>, Set<Path>>>
    adapters: Path?
    phix: Path?

    main:
    QC_MODULE(reads, adapters, phix)

    emit:
    // Individual outputs
    fastq: Channel<Tuple<Map, Set<Path>, Set<Path>>> = QC_MODULE.out.fastq
    fastq_only: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.fastq_only
    txt: Channel<Tuple<Map, Path>> = QC_MODULE.out.txt
    error: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error
    error_fastq: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error_fastq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        QC_MODULE.out.fastq_only,
        QC_MODULE.out.txt,
        QC_MODULE.out.supplemental,
        QC_MODULE.out.error,
        QC_MODULE.out.error_fastq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.versions])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.logs])
}
