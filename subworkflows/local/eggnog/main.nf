//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'eggnog')
options.is_module = params.wf == 'eggnog' ? true : false

downloader_opts = [
    params.skip_diamond ? "-D" : "",
    params.install_hmm ? "-H -d ${hmmer_taxid} " : "",
    params.install_mmseq ? "-M" : "",
    params.install_pfam ? "-P" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()

mapper_opts = [
    "--genepred ${params.genepred}",
    "-m ${params.mode}",
    "${params.eggnog_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { EGGNOG_DOWNLOAD } from '../../../modules/nf-core/eggnog/download/main' addParams( options: options + [ args: "${downloader_opts}"] )
include { EGGNOG_MAPPER } from '../../../modules/nf-core/eggnog/mapper/main' addParams( options:  options +  [ args: "${mapper_opts}"] )

workflow EGGNOG {
    take:
    faa // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    if (params.download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(faa, file("${params.eggnog}/*"))
    }

    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions.first())

    emit:
    hits = EGGNOG_MAPPER.out.hits
    versions = ch_versions
}
