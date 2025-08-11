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
    ch_nf_logs = Channel.empty()

    PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
    PROTEINS = params.proteins ? file(params.proteins) : []

    PROKKA_MODULE(fasta, PROTEINS, PRODIGAL_TF)
    ch_versions = ch_versions.mix(PROKKA_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(PROKKA_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PROKKA_MODULE.out.nf_logs)

    emit:
    annotations = PROKKA_MODULE.out.annotations
    gff = PROKKA_MODULE.out.gff
    gbk = PROKKA_MODULE.out.gbk
    fna = PROKKA_MODULE.out.fna
    faa = PROKKA_MODULE.out.faa
    ffn = PROKKA_MODULE.out.ffn
    sqn = PROKKA_MODULE.out.sqn
    fsa = PROKKA_MODULE.out.fsa
    tbl = PROKKA_MODULE.out.tbl
    txt = PROKKA_MODULE.out.txt
    tsv = PROKKA_MODULE.out.tsv
    blastdb = PROKKA_MODULE.out.blastdb
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
