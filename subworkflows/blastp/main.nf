//
// blastp - Search against protein BLAST databases using protein queries
//
include { BLAST_BLASTP as BLASTP_MODULE } from '../../modules/blast/blastp/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTP {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    QUERY = params.blastp_query ? file(params.blastp_query) : []

    // Run BLASTP
    BLASTP_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTP_MODULE.out.versions)
    ch_logs = ch_logs.mix(BLASTP_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BLASTP_MODULE.out.nf_logs)

    // Merge results
    BLASTP_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'blastp'], tsv] }.set{ ch_merge_blastp }
    CSVTK_CONCAT(ch_merge_blastp, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = BLASTP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
