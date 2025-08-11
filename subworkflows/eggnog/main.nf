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
    ch_nf_logs = Channel.empty()

    DATABASE = params.eggnog_db ? file(params.eggnog_db) : []

    if (params.download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        ch_versions = ch_versions.mix(EGGNOG_DOWNLOAD.out.versions)
        ch_logs = ch_logs.mix(EGGNOG_DOWNLOAD.out.logs)
        ch_nf_logs = ch_nf_logs.mix(EGGNOG_DOWNLOAD.out.nf_logs)
        
        if (params.eggnog_save_as_tarball) {
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db_tarball)
        } else {
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
        }
    } else {
        EGGNOG_MAPPER(faa, DATABASE)
    }
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions.first())
    ch_logs = ch_logs.mix(EGGNOG_MAPPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(EGGNOG_MAPPER.out.nf_logs)

    emit:
    hits = EGGNOG_MAPPER.out.hits
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
