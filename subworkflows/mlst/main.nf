//
// mlst - Automatic MLST calling from assembled contigs
//
include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MLST {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ mlst_db ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    if (params.mlst_db) {
        MLST_MODULE(fasta, file(params.mlst_db))
    } else {
        MLST_MODULE(fasta, db)
    }

    ch_versions = ch_versions.mix(MLST_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(MLST_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MLST_MODULE.out.nf_logs)

    // Merge results
    MLST_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'mlst'], tsv] }.set{ ch_merge_mlst }
    CSVTK_CONCAT(ch_merge_mlst, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = MLST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
