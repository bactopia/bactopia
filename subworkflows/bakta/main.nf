//
// bakta - Rapid annotation of bacterial genomes and plasmids
//
include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN } from '../../modules/bakta/run/main'

workflow BAKTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Set up input files
    DATABASE = params.bakta_db ? file(params.bakta_db) : []
    PROTEINS = params.proteins ? file(params.proteins) : []
    PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
    REPLICONS = params.replicons ? file(params.replicons) : []

    if (params.download_bakta) {
        // Force BAKTA_DOWNLOAD to wait
        BAKTA_DOWNLOAD()

        if (params.bakta_save_as_tarball) {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db_tarball, PROTEINS, PRODIGAL_TF, REPLICONS)
        } else {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db, PROTEINS, PRODIGAL_TF, REPLICONS)
        }
    } else {
        BAKTA_RUN(fasta, DATABASE, PROTEINS, PRODIGAL_TF, REPLICONS)
    }

    ch_versions = ch_versions.mix(BAKTA_RUN.out.versions.first())
    ch_logs = ch_logs.mix(BAKTA_RUN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BAKTA_RUN.out.nf_logs)

    emit:
    annotations = BAKTA_RUN.out.annotations
    embl = BAKTA_RUN.out.embl
    faa = BAKTA_RUN.out.faa
    ffn = BAKTA_RUN.out.ffn
    fna = BAKTA_RUN.out.fna
    gbff = BAKTA_RUN.out.gbff
    gff = BAKTA_RUN.out.gff
    hypotheticals_faa = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv = BAKTA_RUN.out.hypotheticals_tsv
    json = BAKTA_RUN.out.json
    plot = BAKTA_RUN.out.plot
    tsv = BAKTA_RUN.out.tsv
    txt = BAKTA_RUN.out.txt
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
