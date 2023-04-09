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
DATABASE = params.eggnog_db ? file(params.eggnog_db) : []
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
        if (params.eggnog_save_as_tarball) {
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db_tarball)
        } else {
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
        }
    } else {
        EGGNOG_MAPPER(faa, DATABASE)
    }
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions.first())

    emit:
    hits = EGGNOG_MAPPER.out.hits
    versions = ch_versions
}
