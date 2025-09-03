//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'

workflow PROKKA {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]
    proteins
    prodigal_tf

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // Run Prokka
    PROKKA_MODULE(fasta, proteins, prodigal_tf)
    ch_versions = ch_versions.mix(PROKKA_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(PROKKA_MODULE.out.logs)

    emit:
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
    logs = ch_logs
    nf_logs = PROKKA_MODULE.out.nf_begin.mix(
        PROKKA_MODULE.out.nf_err,
        PROKKA_MODULE.out.nf_log,
        PROKKA_MODULE.out.nf_out,
        PROKKA_MODULE.out.nf_run,
        PROKKA_MODULE.out.nf_sh,
        PROKKA_MODULE.out.nf_trace
    )
    versions = ch_versions
}
