//
// prokka - Whole genome annotation of small genomes (bacterial, archeal, viral)
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'prokka')
options.ignore = params.wf == 'pangenome' ? (params.keep_downloads ? [] : [".gz", ".txt", ".tsv"] ) : []
options.logs_subdir = params.wf == 'pangenome' ? "use-prefix" : ""
options.is_module = true
options.args = [
    "--evalue ${params.prokka_evalue}",
    "--coverage ${params.prokka_coverage}",
    "--centre ${params.centre}",
    "${params.prokka_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
PROTEINS = params.proteins ? file(params.proteins) : []
include { PROKKA as PROKKA_MODULE } from '../../../modules/nf-core/prokka/main' addParams( options: options )

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
