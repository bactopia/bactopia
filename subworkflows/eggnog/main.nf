//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER } from '../../modules/eggnog/mapper/main'

workflow EGGNOG {
    take:
    faa // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    if (params.download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        ch_versions = ch_versions.mix(EGGNOG_DOWNLOAD.out.versions)
        ch_logs = ch_logs.mix(EGGNOG_DOWNLOAD.out.logs)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions.first())
    ch_logs = ch_logs.mix(EGGNOG_MAPPER.out.logs)

    emit:
    hits = EGGNOG_MAPPER.out.hits
    DATABASE = params.eggnog_db ? file(params.eggnog_db) : []
    logs = ch_logs
    nf_logs = EGGNOG_DOWNLOAD.out.nf_begin.mix(
        EGGNOG_DOWNLOAD.out.nf_err,
        EGGNOG_DOWNLOAD.out.nf_log,
        EGGNOG_DOWNLOAD.out.nf_out,
        EGGNOG_DOWNLOAD.out.nf_run,
        EGGNOG_DOWNLOAD.out.nf_sh,
        EGGNOG_DOWNLOAD.out.nf_trace
    )
    versions = ch_versions
}
