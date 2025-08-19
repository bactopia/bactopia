//
// blastx - Search against protein BLAST databases using translated nucleotide queries
//
include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTX {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    query // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // Run BLASTX
    BLASTX_MODULE(reads, query)
    ch_versions = ch_versions.mix(BLASTX_MODULE.out.versions)
    ch_logs = ch_logs.mix(BLASTX_MODULE.out.logs)

    // Merge results
    BLASTX_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'blastx'], tsv]}.set{ ch_merge_blastx }
    CSVTK_CONCAT(ch_merge_blastx, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = BLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = BLASTX_MODULE.out.nf_begin.mix(
        BLASTX_MODULE.out.nf_err,
        BLASTX_MODULE.out.nf_log,
        BLASTX_MODULE.out.nf_out,
        BLASTX_MODULE.out.nf_run,
        BLASTX_MODULE.out.nf_sh,
        BLASTX_MODULE.out.nf_trace,
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
