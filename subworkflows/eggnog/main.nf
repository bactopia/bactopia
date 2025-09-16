//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER } from '../../modules/eggnog/mapper/main'

workflow EGGNOG {
    take:
    faa // channel: [ val(meta), [ fasta ] ]
    database
    download_eggnog

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    if (download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(faa, database)
    }
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)
    ch_logs = ch_logs.mix(EGGNOG_MAPPER.out.logs)

    emit:
    hits = EGGNOG_MAPPER.out.hits
    logs = ch_logs
    nf_logs = EGGNOG_MAPPER.out.nf_begin.mix(
        EGGNOG_MAPPER.out.nf_err,
        EGGNOG_MAPPER.out.nf_log,
        EGGNOG_MAPPER.out.nf_out,
        EGGNOG_MAPPER.out.nf_run,
        EGGNOG_MAPPER.out.nf_sh,
        EGGNOG_MAPPER.out.nf_trace
    )
    versions = ch_versions
}
