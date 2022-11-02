//
// kleborate - Screening Klebsiella genome assemblies for MLST, sub-species, and other Klebsiella related genes of interest
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'kleborate')
options.is_module = params.wf == 'kleborate' ? true : false
options.args = [
    params.skip_resistance ? "" : "--resistance",
    params.skip_kaptive ? "" : "--kaptive",
    params.force_index ? "--force_index" : "",
    "--min_identity ${params.min_identity}",
    "--min_coverage ${params.kleborate_min_coverage}",
    "--min_spurious_identity ${params.min_spurious_identity}",
    "--min_spurious_coverage ${params.min_spurious_coverage}",
    "--min_kaptive_confidence ${params.min_kaptive_confidence}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { KLEBORATE as KLEBORATE_MODULE } from '../../../modules/nf-core/kleborate/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'kleborate'] )
}

workflow KLEBORATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_kleborate = Channel.empty()

    KLEBORATE_MODULE(fasta)
    ch_versions = ch_versions.mix(KLEBORATE_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        KLEBORATE_MODULE.out.txt.collect{meta, txt -> txt}.map{ txt -> [[id:'kleborate'], txt]}.set{ ch_merge_kleborate }
        CSVTK_CONCAT(ch_merge_kleborate, 'tsv', 'tsv')
        ch_merged_kleborate = ch_merged_kleborate.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = KLEBORATE_MODULE.out.txt
    merged_tsv = ch_merged_kleborate
    versions = ch_versions
}
