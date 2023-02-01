//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
prokka_args = [
    params.compliant ? "--compliant" : "",
    "--evalue ${params.prokka_evalue}",
    "--coverage ${params.prokka_coverage}",
    "${params.prokka_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
PROTEINS = params.proteins ? file(params.proteins) : []

include { PROKKA as PROKKA_MODULE } from '../../../modules/nf-core/prokka/main' addParams( options: [ args: "${prokka_args}", is_module: true] )
include { PROKKA_MAIN as USE_PROKKA } from '../../../modules/nf-core/prokka/main' addParams( options: [ args: "${prokka_args}", is_module: true] )

workflow PROKKA {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    PROKKA_MODULE(fasta, PROTEINS, PRODIGAL_TF)
    ch_versions = ch_versions.mix(PROKKA_MODULE.out.versions.first())

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
    versions = ch_versions // channel: [ versions.yml ]
}

workflow PROKKA_MAIN {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    USE_PROKKA(fasta)
    ch_versions = ch_versions.mix(USE_PROKKA.out.versions.first())

    emit:
    annotations = USE_PROKKA.out.annotations
    gff = USE_PROKKA.out.gff
    gbk = USE_PROKKA.out.gbk
    fna = USE_PROKKA.out.fna
    faa = USE_PROKKA.out.faa
    ffn = USE_PROKKA.out.ffn
    sqn = USE_PROKKA.out.sqn
    fsa = USE_PROKKA.out.fsa
    tbl = USE_PROKKA.out.tbl
    txt = USE_PROKKA.out.txt
    tsv = USE_PROKKA.out.tsv
    blastdb = PROKKA_MODULE.out.blastdb
    versions = ch_versions // channel: [ versions.yml ]
}
