//
// midas - Estimate species abundances from FASTQ files
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'midas')
options.is_module = params.wf == 'midas' ? true : false
options.args = [
    "--word_size ${params.midas_word_size}",
    "--aln_cov ${params.midas_aln_cov}",
    params.midas_debug ? "" : "--remove_temp",
    "${params.midas_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE = params.midas_db ? file(params.midas_db) : []

include { MIDAS_SPECIES } from '../../../modules/nf-core/midas/species/main' addParams( options: options )

workflow MIDAS {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MIDAS_SPECIES(reads, DATABASE)
    ch_versions = ch_versions.mix(MIDAS_SPECIES.out.versions)

    emit:
    tsv = MIDAS_SPECIES.out.tsv
    abundances = MIDAS_SPECIES.out.abundances
    versions = ch_versions // channel: [ versions.yml ]
}
