//
// clermontyping - in silico phylotyping of Escherichia genus
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'clermontyping')
options.args = [
    "--threshold ${params.clermon_threshold}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../../modules/nf-core/clermontyping/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'clermontyping-concat', process_name: params.merge_folder] )

workflow CLERMONTYPING {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    CLERMONTYPING_MODULE(fasta)
    ch_versions = ch_versions.mix(CLERMONTYPING_MODULE.out.versions.first())

    // Merge results
    CLERMONTYPING_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'clermontyping'], tsv]}.set{ ch_merge_clermontyping }
    CSVTK_CONCAT(ch_merge_clermontyping, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = CLERMONTYPING_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
