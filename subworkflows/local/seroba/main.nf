//
// seroba - Serotyping of Streptococcus pneumoniae from sequence reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'seroba')
options.is_module = params.wf == 'seroba' ? true : false
options.args = [
    "--coverage ${params.seroba_coverage}",
    params.seroba_noclean ? "--noclean" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { SEROBA_RUN } from '../../../modules/nf-core/modules/seroba/run/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'seroba'] )
}

workflow SEROBA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_seroba = Channel.empty()

    SEROBA_RUN(fasta)
    ch_versions = ch_versions.mix(SEROBA_RUN.out.versions.first())

    if (params.is_subworkflow) {
        SEROBA_RUN.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'seroba'], tsv]}.set{ ch_merge_seroba }
        CSVTK_CONCAT(ch_merge_seroba, 'tsv', 'tsv')
        ch_merged_seroba = ch_merged_seroba.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SEROBA_RUN.out.tsv
    merged_tsv = ch_merged_seroba
    versions = ch_versions
}
