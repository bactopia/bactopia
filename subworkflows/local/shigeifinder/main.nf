//
// shigeifinder - Shigella and EIEC serotyping from assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigeifinder')
options.is_module = params.wf == 'shigeifinder' ? true : false

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../../modules/nf-core/shigeifinder/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'shigeifinder'] )
}

workflow SHIGEIFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_shigeifinder= Channel.empty()

    SHIGEIFINDER_MODULE(fasta)
    ch_versions = ch_versions.mix(SHIGEIFINDER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SHIGEIFINDER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'shigeifinder'], tsv]}.set{ ch_merge_shigeifinder }
        CSVTK_CONCAT(ch_merge_shigeifinder, 'tsv', 'tsv')
        ch_merged_shigeifinder = ch_merged_shigeifinder.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SHIGEIFINDER_MODULE.out.tsv
    merged_tsv = ch_merged_shigeifinder
    versions = ch_versions
}
