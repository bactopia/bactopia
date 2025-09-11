//
// shigeifinder - Shigella and EIEC serotyping from assemblies
//
include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../modules/shigeifinder/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SHIGEIFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SHIGEIFINDER_MODULE(fasta)
    ch_versions = ch_versions.mix(SHIGEIFINDER_MODULE.out.versions)

    // Merge results
    SHIGEIFINDER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'shigeifinder'], tsv]}.set{ ch_merge_shigeifinder }
    CSVTK_CONCAT(ch_merge_shigeifinder, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SHIGEIFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = SHIGEIFINDER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SHIGEIFINDER_MODULE.out.nf_begin.mix(
        SHIGEIFINDER_MODULE.out.nf_err,
        SHIGEIFINDER_MODULE.out.nf_log,
        SHIGEIFINDER_MODULE.out.nf_out,
        SHIGEIFINDER_MODULE.out.nf_run,
        SHIGEIFINDER_MODULE.out.nf_sh,
        SHIGEIFINDER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions
}
