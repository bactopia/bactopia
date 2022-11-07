//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pasty')
options.is_module = params.wf == 'pasty' ? true : false
options.args = [
    "--min_pident ${params.pasty_min_pident}",
    "--min_coverage ${params.pasty_min_coverage}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PASTY as PASTY_MODULE } from '../../../modules/nf-core/pasty/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'pasty'] )
}

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_pasty = Channel.empty()

    PASTY_MODULE(fasta)
    ch_versions = ch_versions.mix(PASTY_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        PASTY_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'pasty'], tsv]}.set{ ch_merge_pasty }
        CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')
        ch_merged_pasty = ch_merged_pasty.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = ch_merged_pasty
    blast = PASTY_MODULE.out.blast
    versions = ch_versions
}
