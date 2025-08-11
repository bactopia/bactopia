//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    BLASTDB = params.emmtyper_blastdb ? file(params.emmtyper_blastdb) : []

    EMMTYPER_MODULE(fasta, BLASTDB)
    ch_versions = ch_versions.mix(EMMTYPER_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(EMMTYPER_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(EMMTYPER_MODULE.out.nf_logs)

    // Merge results
    EMMTYPER_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'emmtyper'], tsv] }.set{ ch_merge_emmtyper }
    CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
