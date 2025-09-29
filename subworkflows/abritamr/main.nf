//
// abritamr - A NATA accredited tool for reporting the presence of antimicrobial resistance genes
//
include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ABRITAMR {
    take:
    fasta // channel: [ val(meta), [ reads ] ]

    main:
    ABRITAMR_RUN ( fasta )

    // Merge results
    ABRITAMR_RUN.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:'abritamr'], summary]}.set{ ch_merge_summary }
    CSVTK_CONCAT(ch_merge_summary, 'tsv', 'tsv')

    emit:
    // Individual outputs
    summary_tsv = ABRITAMR_RUN.out.summary
    merged_summary_tsv = CSVTK_CONCAT.out.csv
    matches_tsv = ABRITAMR_RUN.out.matches
    partials_tsv = ABRITAMR_RUN.out.partials
    virulence_tsv = ABRITAMR_RUN.out.virulence
    amrfinder_tsv = ABRITAMR_RUN.out.amrfinder

    // Generic aggregate outputs
    results = ABRITAMR_RUN.out.summary.mix(
        CSVTK_CONCAT.out.csv,
        ABRITAMR_RUN.out.matches,
        ABRITAMR_RUN.out.partials,
        ABRITAMR_RUN.out.virulence,
        ABRITAMR_RUN.out.amrfinder
    )
    logs = ABRITAMR_RUN.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = ABRITAMR_RUN.out.nf_begin.mix(
        ABRITAMR_RUN.out.nf_err,
        ABRITAMR_RUN.out.nf_log,
        ABRITAMR_RUN.out.nf_out,
        ABRITAMR_RUN.out.nf_run,
        ABRITAMR_RUN.out.nf_sh,
        ABRITAMR_RUN.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ABRITAMR_RUN.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
