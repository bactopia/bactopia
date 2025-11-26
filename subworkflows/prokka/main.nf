//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
nextflow.preview.types = true

include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PROKKA {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ assemblies ] ]
    proteins: Channel<Tuple<Map, Path>>
    prodigal_tf: Channel<Tuple<Map, Path>>

    main:
    PROKKA_MODULE(fasta, proteins, prodigal_tf)

    emit:
    // Individual outputs
    annotations: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.annotations
    blastdb: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.blastdb
    faa: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.faa
    ffn: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.ffn
    fna: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.fna
    fsa: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.fsa
    gbk: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.gbk
    gff: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.gff
    tsv: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.tsv
    txt: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.txt
    sqn: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.sqn
    tbl: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.tbl

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PROKKA_MODULE.out.blastdb,
        PROKKA_MODULE.out.faa,
        PROKKA_MODULE.out.ffn,
        PROKKA_MODULE.out.fna,
        PROKKA_MODULE.out.fsa,
        PROKKA_MODULE.out.gbk,
        PROKKA_MODULE.out.gff,
        PROKKA_MODULE.out.tsv,
        PROKKA_MODULE.out.txt,
        PROKKA_MODULE.out.sqn,
        PROKKA_MODULE.out.tbl
    ])
    logs: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.versions
}
