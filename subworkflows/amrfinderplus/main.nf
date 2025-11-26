//
// amrfinderplus - Identify antimicrobial resistance in genes or proteins
//
nextflow.preview.types = true

include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT      } from '../../modules/csvtk/concat/main'
include { flattenPaths      } from 'plugin/nf-bactopia'
include { gather            } from 'plugin/nf-bactopia'

workflow AMRFINDERPLUS {
    take:
    fasta: Channel<Tuple<Map, Path>>
    db: Channel<Tuple<Map, Path>>

    main:
    AMRFINDERPLUS_RUN(fasta, db)
    CSVTK_CONCAT(gather(AMRFINDERPLUS_RUN.out.report, 'amrfinderplus'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = AMRFINDERPLUS_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    mutation_report: Channel<Tuple<Map, Path>> = AMRFINDERPLUS_RUN.out.mutation_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.report,
        AMRFINDERPLUS_RUN.out.mutation_report,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
