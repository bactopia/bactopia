//
// abritamr - A NATA accredited tool for reporting the presence of antimicrobial resistance genes
//
nextflow.preview.types = true

include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow ABRITAMR {
    take:
    fasta: Channel<Tuple<Map, Path>>

    main:
    ABRITAMR_RUN(fasta)
    CSVTK_CONCAT(gather(ABRITAMR_RUN.out.summary, 'abritamr'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    summary_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.summary
    merged_summary_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    matches_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.matches
    partials_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.partials
    virulence_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.virulence
    amrfinder_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.amrfinder

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.summary,
        ABRITAMR_RUN.out.matches,
        ABRITAMR_RUN.out.partials,
        ABRITAMR_RUN.out.virulence,
        ABRITAMR_RUN.out.amrfinder,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
