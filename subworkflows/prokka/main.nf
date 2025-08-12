//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'

workflow PROKKA {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    PROKKA_MODULE(fasta, PROTEINS, PRODIGAL_TF)
    ch_versions = ch_versions.mix(PROKKA_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(PROKKA_MODULE.out.logs)

    emit:
    tsv = PROKKA_MODULE.out.tsv
    txt = PROKKA_MODULE.out.txt
    gbk = PROKKA_MODULE.out.gbk
    PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
    PROTEINS = params.proteins ? file(params.proteins) : []
    annotations = PROKKA_MODULE.out.annotations
    blastdb = PROKKA_MODULE.out.blastdb
    faa = PROKKA_MODULE.out.faa
    ffn = PROKKA_MODULE.out.ffn
    fna = PROKKA_MODULE.out.fna
    fsa = PROKKA_MODULE.out.fsa
    gff = PROKKA_MODULE.out.gff
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
