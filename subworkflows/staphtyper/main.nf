//
// staphtyper - Determine the agr, spa and SCCMEC types for S. aureus assemblies
//
include { AGRVATE } from '../../modules/agrvate/main'
include { SPATYPER } from '../../modules/spatyper/main'
include { SCCMEC } from '../../modules/sccmec/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_AGRVATE } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SPATYPER } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SCCMEC } from '../../modules/csvtk/concat/main'

workflow STAPHTYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    REPEATS = params.repeats ? file(params.repeats, checkIfExists: true) : []
    REPEAT_ORDER = params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []
    
    // agrvate - agr locus type and agr operon variants
    // spatyper - spa typing
    // sccmec - SCCmec type based on targets and full cassettes
    AGRVATE(fasta)
    SPATYPER(fasta, REPEATS, REPEAT_ORDER)
    SCCMEC(fasta)

    // gather versions
    ch_versions = ch_versions.mix(AGRVATE.out.versions)
    ch_versions = ch_versions.mix(SPATYPER.out.versions)
    ch_versions = ch_versions.mix(SCCMEC.out.versions)

    // Merge AgrVATE
    AGRVATE.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:'agrvate'], summary]}.set{ ch_merge_agrvate }
    CSVTK_CONCAT_AGRVATE(ch_merge_agrvate, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_AGRVATE.out.versions)

    // Merge spaTyper
    SPATYPER.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'spatyper'], tsv]}.set{ ch_merge_spatyper }
    CSVTK_CONCAT_SPATYPER(ch_merge_spatyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SPATYPER.out.versions)

    // Merge sccmec
    SCCMEC.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'sccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT_SCCMEC(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SCCMEC.out.versions)

    emit:
    agrvate_tsv = AGRVATE.out.summary
    agrvate_merged_tsv = CSVTK_CONCAT_AGRVATE.out.csv
    spatyper_tsv = SPATYPER.out.tsv
    spatyper_merged_tsv = CSVTK_CONCAT_SPATYPER.out.csv
    sccmec_tsv = SCCMEC.out.tsv
    sccmec_merged_tsv = CSVTK_CONCAT_SCCMEC.out.csv
    logs = AGRVATE.out.logs.mix(
        SPATYPER.out.logs,
        SCCMEC.out.logs,
        CSVTK_CONCAT_AGRVATE.out.logs,
        CSVTK_CONCAT_SPATYPER.out.logs,
        CSVTK_CONCAT_SCCMEC.out.logs
    )
    nf_logs = AGRVATE.out.nf_begin.mix(
        AGRVATE.out.nf_err,
        AGRVATE.out.nf_log,
        AGRVATE.out.nf_out,
        AGRVATE.out.nf_run,
        AGRVATE.out.nf_sh,
        AGRVATE.out.nf_trace,
        SPATYPER.out.nf_begin,
        SPATYPER.out.nf_err,
        SPATYPER.out.nf_log,
        SPATYPER.out.nf_out,
        SPATYPER.out.nf_run,
        SPATYPER.out.nf_sh,
        SPATYPER.out.nf_trace,
        SCCMEC.out.nf_begin,
        SCCMEC.out.nf_err,
        SCCMEC.out.nf_log,
        SCCMEC.out.nf_out,
        SCCMEC.out.nf_run,
        SCCMEC.out.nf_sh,
        SCCMEC.out.nf_trace,
        CSVTK_CONCAT_AGRVATE.out.nf_begin,
        CSVTK_CONCAT_AGRVATE.out.nf_err,
        CSVTK_CONCAT_AGRVATE.out.nf_log,
        CSVTK_CONCAT_AGRVATE.out.nf_out,
        CSVTK_CONCAT_AGRVATE.out.nf_run,
        CSVTK_CONCAT_AGRVATE.out.nf_sh,
        CSVTK_CONCAT_AGRVATE.out.nf_trace,
        CSVTK_CONCAT_SPATYPER.out.nf_begin,
        CSVTK_CONCAT_SPATYPER.out.nf_err,
        CSVTK_CONCAT_SPATYPER.out.nf_log,
        CSVTK_CONCAT_SPATYPER.out.nf_out,
        CSVTK_CONCAT_SPATYPER.out.nf_run,
        CSVTK_CONCAT_SPATYPER.out.nf_sh,
        CSVTK_CONCAT_SPATYPER.out.nf_trace,
        CSVTK_CONCAT_SCCMEC.out.nf_begin,
        CSVTK_CONCAT_SCCMEC.out.nf_err,
        CSVTK_CONCAT_SCCMEC.out.nf_log,
        CSVTK_CONCAT_SCCMEC.out.nf_out,
        CSVTK_CONCAT_SCCMEC.out.nf_run,
        CSVTK_CONCAT_SCCMEC.out.nf_sh,
        CSVTK_CONCAT_SCCMEC.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
