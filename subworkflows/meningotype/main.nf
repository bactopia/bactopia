//
// meningotype - Serotyping of Neisseria meningitidis
//
include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    MENINGOTYPE_MODULE(fasta)
    ch_versions = ch_versions.mix(MENINGOTYPE_MODULE.out.versions)
    ch_logs = ch_logs.mix(MENINGOTYPE_MODULE.out.logs)

    // Merge results
    MENINGOTYPE_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'meningotype'], tsv]}.set{ ch_merge_meningotype }
    CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        MENINGOTYPE_MODULE.out.nf_begin,
        MENINGOTYPE_MODULE.out.nf_err,
        MENINGOTYPE_MODULE.out.nf_log,
        MENINGOTYPE_MODULE.out.nf_out,
        MENINGOTYPE_MODULE.out.nf_run,
        MENINGOTYPE_MODULE.out.nf_sh,
        MENINGOTYPE_MODULE.out.nf_trace
    )
    versions = ch_versions
}
