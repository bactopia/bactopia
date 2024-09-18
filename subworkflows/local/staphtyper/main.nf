//
// staphtyper - Determine the agr, spa and SCCMEC types for S. aureus assemblies
//
REPEATS = params.repeats ? file(params.repeats, checkIfExists: true) : []
REPEAT_ORDER = params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []

argvate_args = params.typing_only ? '--typing_only' : ''
spatyper_args = params.do_enrich ? '--do_enrich' : ''

include { AGRVATE } from '../../../modules/nf-core/agrvate/main' addParams( options: [args: "${argvate_args}", is_subworkflow: true] )
include { SPATYPER } from '../../../modules/nf-core/spatyper/main' addParams( options: [args: "${spatyper_args}", is_subworkflow: true] )
include { SCCMEC } from '../../../modules/nf-core/sccmec/main' addParams( options: [is_subworkflow: true] )

include { CSVTK_CONCAT as CSVTK_CONCAT_AGRVATE } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [args: '-C "$"', logs_subdir: 'agrvate-concat', process_name: params.merge_folder])
include { CSVTK_CONCAT as CSVTK_CONCAT_SPATYPER } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'spatyper-concat', process_name: params.merge_folder])
include { CSVTK_CONCAT as CSVTK_CONCAT_SCCMEC } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'sccmec-concat', process_name: params.merge_folder])

workflow STAPHTYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    // agrvate - agr locus type and agr operon variants
    // spatyper - spa typing
    // sccmec - SCCmec type based on targets and full cassettes
    AGRVATE(fasta)
    SPATYPER(fasta, REPEATS, REPEAT_ORDER)
    SCCMEC(fasta)

    // gather versions
    ch_versions = ch_versions.mix(AGRVATE.out.versions.first())
    ch_versions = ch_versions.mix(SPATYPER.out.versions.first())
    ch_versions = ch_versions.mix(SCCMEC.out.versions.first())

    // Merge AgrVATE
    AGRVATE.out.summary.collect{meta, summary -> summary}.map{ summary -> [[id:'agrvate'], summary]}.set{ ch_merge_agrvate }
    CSVTK_CONCAT_AGRVATE(ch_merge_agrvate, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_AGRVATE.out.versions)

    // Merge spaTyper
    SPATYPER.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'spatyper'], tsv]}.set{ ch_merge_spatyper }
    CSVTK_CONCAT_SPATYPER(ch_merge_spatyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SPATYPER.out.versions)

    // Merge sccmec
    SCCMEC.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'sccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT_SCCMEC(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SCCMEC.out.versions)

    emit:
    agrvate_tsv = AGRVATE.out.summary
    agrvate_merged_tsv = CSVTK_CONCAT_AGRVATE.out.csv
    spatyper_tsv = SPATYPER.out.tsv
    spatyper_merged_tsv = CSVTK_CONCAT_SPATYPER.out.csv
    sccmec_tsv = SCCMEC.out.tsv
    sccmec_merged_tsv = CSVTK_CONCAT_SCCMEC.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
