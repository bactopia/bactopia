//
// abritamr - A NATA accredited tool for reporting the presence of antimicrobial resistance genes
//
nextflow.preview.types = true

include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ABRITAMR {
    take:
    fasta: Channel<Tuple<Map, Path>>

    main:
    ABRITAMR_RUN(fasta)

    // Merge results
    ch_merge_summary = ABRITAMR_RUN.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:'abritamr'], summary]}
    CSVTK_CONCAT(ch_merge_summary, 'tsv', 'tsv')

    emit:
    // Individual outputs
    summary_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.summary
    merged_summary_tsv: Channel<Tuple<Map, Path>> = channel.of(CSVTK_CONCAT.out.csv)
    matches_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.matches
    partials_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.partials
    virulence_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.virulence
    amrfinder_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.amrfinder

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.summary.mix(
        ABRITAMR_RUN.out.matches,
        ABRITAMR_RUN.out.partials,
        ABRITAMR_RUN.out.virulence,
        ABRITAMR_RUN.out.amrfinder,
        channel.of(CSVTK_CONCAT.out.csv)
    )
    logs: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.logs.mix(
        channel.of(CSVTK_CONCAT.out.logs)
    )
    nf_logs: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.nf_begin.mix(
        ABRITAMR_RUN.out.nf_err,
        ABRITAMR_RUN.out.nf_log,
        ABRITAMR_RUN.out.nf_out,
        ABRITAMR_RUN.out.nf_run,
        ABRITAMR_RUN.out.nf_sh,
        ABRITAMR_RUN.out.nf_trace,
        channel.of(CSVTK_CONCAT.out.nf_begin),
        channel.of(CSVTK_CONCAT.out.nf_err),
        channel.of(CSVTK_CONCAT.out.nf_log),
        channel.of(CSVTK_CONCAT.out.nf_out),
        channel.of(CSVTK_CONCAT.out.nf_run),
        channel.of(CSVTK_CONCAT.out.nf_sh),
        channel.of(CSVTK_CONCAT.out.nf_trace)
    )
    versions: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.versions.mix(
        channel.of(CSVTK_CONCAT.out.versions)
    )
}
