//
// kleborate - Screening Klebsiella genome assemblies for MLST, sub-species, and other Klebsiella related genes of interest
//
include { KLEBORATE as KLEBORATE_MODULE } from '../../modules/kleborate/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow KLEBORATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    KLEBORATE_MODULE(fasta)
    ch_versions = ch_versions.mix(KLEBORATE_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(KLEBORATE_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KLEBORATE_MODULE.out.nf_logs)

    // Merge results
    KLEBORATE_MODULE.out.txt.collect{ _meta, txt -> txt }.map{ txt -> [[id:'kleborate'], txt] }.set{ ch_merge_kleborate }
    CSVTK_CONCAT(ch_merge_kleborate, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = KLEBORATE_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
