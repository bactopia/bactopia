//
// bakta - Rapid annotation of bacterial genomes and plasmids
//
nextflow.preview.types = true

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

    emit:
    // Individual outputs
    tsv = BAKTA_RUN.out.tsv
    txt = BAKTA_RUN.out.txt
    embl = BAKTA_RUN.out.embl
    faa = BAKTA_RUN.out.faa
    ffn = BAKTA_RUN.out.ffn
    fna = BAKTA_RUN.out.fna
    gbff = BAKTA_RUN.out.gbff
    gff = BAKTA_RUN.out.gff
    hypotheticals_faa = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv = BAKTA_RUN.out.hypotheticals_tsv
    blastdb = BAKTA_RUN.out.blastdb
    annotations = BAKTA_RUN.out.annotations

    // Generic aggregate outputs
    results = BAKTA_RUN.out.tsv.mix(
        BAKTA_RUN.out.txt,
        BAKTA_RUN.out.embl,
        BAKTA_RUN.out.faa,
        BAKTA_RUN.out.ffn,
        BAKTA_RUN.out.fna,
        BAKTA_RUN.out.gbff,
        BAKTA_RUN.out.gff,
        BAKTA_RUN.out.hypotheticals_faa,
        BAKTA_RUN.out.hypotheticals_tsv,
        BAKTA_RUN.out.blastdb
    )
    logs = BAKTA_RUN.out.logs
    nf_logs = BAKTA_RUN.out.nf_begin.mix(
        BAKTA_RUN.out.nf_err,
        BAKTA_RUN.out.nf_log,
        BAKTA_RUN.out.nf_out,
        BAKTA_RUN.out.nf_run,
        BAKTA_RUN.out.nf_sh,
        BAKTA_RUN.out.nf_trace
    )
    versions = BAKTA_RUN.out.versions
}
