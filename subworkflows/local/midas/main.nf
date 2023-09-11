//
// midas - Estimate species abundances from FASTQ files
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'midas')
options.args = [
    "--word_size ${params.midas_word_size}",
    "--aln_cov ${params.midas_aln_cov}",
    params.midas_debug ? "" : "--remove_temp",
    "${params.midas_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE = params.midas_db ? file(params.midas_db) : []

include { MIDAS_SPECIES } from '../../../modules/nf-core/midas/species/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'midas-concat', process_name: params.merge_folder] )

workflow MIDAS {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MIDAS_SPECIES(reads, DATABASE)
    ch_versions = ch_versions.mix(MIDAS_SPECIES.out.versions)

    // Merge results
    MIDAS_SPECIES.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'midas'], tsv]}.set{ ch_merge_midas }
    CSVTK_CONCAT(ch_merge_midas, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = MIDAS_SPECIES.out.tsv
    abundances = MIDAS_SPECIES.out.abundances
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
