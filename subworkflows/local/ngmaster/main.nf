//
// ngmaster - Multi-antigen sequence typing for Neisseria gonorrhoeae
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ngmaster')
options.is_module = params.wf == 'ngmaster' ? true : false
options.args = [
    params.csv ? "--csv" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { NGMASTER as NGMASTER_MODULE } from '../../../modules/nf-core/ngmaster/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'ngmaster'] )
}

workflow NGMASTER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_ngmaster = Channel.empty()

    NGMASTER_MODULE(fasta)
    ch_versions = ch_versions.mix(NGMASTER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        NGMASTER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ngmaster'], tsv]}.set{ ch_merge_ngmaster }
        CSVTK_CONCAT(ch_merge_ngmaster, 'tsv', 'tsv')
        ch_merged_ngmaster = ch_merged_ngmaster.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = NGMASTER_MODULE.out.tsv
    merged_tsv = ch_merged_ngmaster
    versions = ch_versions
}
