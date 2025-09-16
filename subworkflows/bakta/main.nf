//
// bakta - Rapid annotation of bacterial genomes and plasmids
//
include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN } from '../../modules/bakta/run/main'

workflow BAKTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    database
    download_bakta
    save_as_tarball
    proteins
    prodigal_tf
    replicons


    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    // Set up input files


    if (download_bakta) {
        // Force BAKTA_DOWNLOAD to wait
        BAKTA_DOWNLOAD()

        if (save_as_tarball) {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db_tarball, proteins, prodigal_tf, replicons)
        } else {
            BAKTA_RUN(fasta, BAKTA_DOWNLOAD.out.db, proteins, prodigal_tf, replicons)
        }
    } else {
        BAKTA_RUN(fasta, database, proteins, prodigal_tf, replicons)
    }

    ch_versions = ch_versions.mix(BAKTA_RUN.out.versions)
    ch_logs = ch_logs.mix(BAKTA_RUN.out.logs)

    emit:
    tsv = BAKTA_RUN.out.tsv
    txt = BAKTA_RUN.out.txt
    annotations = BAKTA_RUN.out.annotations
    embl = BAKTA_RUN.out.embl
    faa = BAKTA_RUN.out.faa
    ffn = BAKTA_RUN.out.ffn
    fna = BAKTA_RUN.out.fna
    gbff = BAKTA_RUN.out.gbff
    gff = BAKTA_RUN.out.gff
    hypotheticals_faa = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv = BAKTA_RUN.out.hypotheticals_tsv
    logs = ch_logs
    nf_logs = BAKTA_RUN.out.nf_begin.mix(
        BAKTA_RUN.out.nf_err,
        BAKTA_RUN.out.nf_log,
        BAKTA_RUN.out.nf_out,
        BAKTA_RUN.out.nf_run,
        BAKTA_RUN.out.nf_sh,
        BAKTA_RUN.out.nf_trace
    )
    versions = ch_versions
}
