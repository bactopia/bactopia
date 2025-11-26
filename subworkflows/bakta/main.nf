//
// bakta - Rapid annotation of bacterial genomes and plasmids
//
nextflow.preview.types = true

include { BAKTA_DOWNLOAD } from '../../modules/bakta/download/main'
include { BAKTA_RUN      } from '../../modules/bakta/run/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow BAKTA {
    take:
    fasta: Channel<Tuple<Map, Path>>
    database: Channel<Tuple<Map, Path>>
    download_bakta: Channel<Tuple<Map, Path>>
    save_as_tarball: Channel<Tuple<Map, Path>>
    proteins: Channel<Tuple<Map, Path>>
    prodigal_tf: Channel<Tuple<Map, Path>>
    replicons: Channel<Tuple<Map, Path>>

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
    tsv: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.tsv
    txt: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.txt
    embl: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.embl
    faa: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.faa
    ffn: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.ffn
    fna: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.fna
    gbff: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.gbff
    gff: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.gff
    hypotheticals_faa: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.hypotheticals_faa
    hypotheticals_tsv: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.hypotheticals_tsv
    blastdb: Channel<Tuple<Map, Set<Path>>> = BAKTA_RUN.out.blastdb
    annotations: Channel<Tuple<Map, Path>> = BAKTA_RUN.out.annotations

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BAKTA_RUN.out.tsv,
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
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([BAKTA_RUN.out.versions])
}
