//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
downloader_opts = [
    params.skip_diamond ? "-D" : "",
    params.install_hmm ? "-H -d ${hmmer_taxid} " : "",
    params.install_mmseq ? "-M" : "",
    params.install_pfam ? "-P" : "",=
].join(' ').replaceAll("\\s{2,}", " ").trim()

mapper_opts = [
    "--genepred ${params.genepred}?",
    "-m ${params.mode}",
    "${params.eggnog_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { EGGNOG_DOWNLOAD } from '../../../modules/nf-core/modules/eggnog/download/main' addParams( options: [ args: "${downloader_opts}"] )
include { EGGNOG_MAPPER } from '../../../modules/nf-core/modules/eggnog/mapper/main' addParams( options: [ args: "${mapper_opts}"] )

workflow EGGNOG {
    take:
    faa // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    EGGNOG_DOWNLOAD()
    EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions.first())

    emit:
    tsv = EGGNOG_MAPPER.out.tsv
    versions = ch_versions
}
