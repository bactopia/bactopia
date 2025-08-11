//
// blastx - Search against protein BLAST databases using translated nucleotide queries
//
include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTX {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    QUERY = params.blastx_query ? file(params.blastx_query) : []

    // Run BLASTX
    BLASTX_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTX_MODULE.out.versions)
    ch_logs = ch_logs.mix(BLASTX_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BLASTX_MODULE.out.nf_logs)

    // Merge results
    BLASTX_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'blastx'], tsv] }.set{ ch_merge_blastx }
    CSVTK_CONCAT(ch_merge_blastx, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = BLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
