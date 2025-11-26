//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
nextflow.preview.types = true

include { ABRICATE_RUN     } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow ABRICATE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    ABRICATE_RUN(fasta)
    ABRICATE_SUMMARY(gather(ABRICATE_RUN.out.report, 'abricate'))

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = ABRICATE_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Path>> = ABRICATE_SUMMARY.out.report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = ABRICATE_RUN.out.report.mix(
        ABRICATE_SUMMARY.out.report
    )
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.logs,
        ABRICATE_SUMMARY.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.nf_logs,
        ABRICATE_SUMMARY.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = ABRICATE_RUN.out.versions.mix(
        ABRICATE_SUMMARY.out.versions
    )
}
