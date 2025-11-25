//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
nextflow.preview.types = true

include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'

workflow PROKKA {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]
    proteins
    prodigal_tf

    main:
    PROKKA_MODULE(fasta, proteins, prodigal_tf)

    emit:
    // Individual outputs
    annotations = PROKKA_MODULE.out.annotations
    blastdb = PROKKA_MODULE.out.blastdb
    faa = PROKKA_MODULE.out.faa
    ffn = PROKKA_MODULE.out.ffn
    fna = PROKKA_MODULE.out.fna
    fsa = PROKKA_MODULE.out.fsa
    gbk = PROKKA_MODULE.out.gbk
    gff = PROKKA_MODULE.out.gff
    tsv = PROKKA_MODULE.out.tsv
    txt = PROKKA_MODULE.out.txt
    sqn = PROKKA_MODULE.out.sqn
    tbl = PROKKA_MODULE.out.tbl

    // Generic aggregate outputs
    results = PROKKA_MODULE.out.blastdb.mix(
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
    )
    logs = PROKKA_MODULE.out.logs
    nf_logs = PROKKA_MODULE.out.nf_begin.mix(
        PROKKA_MODULE.out.nf_err,
        PROKKA_MODULE.out.nf_log,
        PROKKA_MODULE.out.nf_out,
        PROKKA_MODULE.out.nf_run,
        PROKKA_MODULE.out.nf_sh,
        PROKKA_MODULE.out.nf_trace
    )
    versions = PROKKA_MODULE.out.versions
}
