//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
PROTEINS = params.proteins ? file(params.proteins) : []
prokka_args = [
    params.compliant ? "--compliant" : "",
    "--evalue ${params.prokka_evalue}",
    "--coverage ${params.prokka_coverage}",
    "${params.prokka_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PROKKA as PROKKA_MODULE } from '../../../modules/nf-core/prokka/main' addParams( options: [ args: "${prokka_args}", is_module: true] )

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
    versions = ch_versions // channel: [ versions.yml ]
}
