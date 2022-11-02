//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'emmtyper')
options.is_module = params.wf == 'emmtyper' ? true : false
options.args = [
    "--workflow ${params.emmtyper_wf}",
    "--cluster-distance ${params.cluster_distance}",
    "--percent-identity ${params.percid}",
    "--culling-limit ${params.culling_limit}",
    "--mismatch ${params.mismatch}",
    "--align-diff ${params.align_diff}",
    "--gap ${params.gap}",
    "--min-perfect ${params.min_perfect}",
    "--min-good ${params.min_good}",
    "--max-size ${params.max_size}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { EMMTYPER as EMMTYPER_MODULE } from '../../../modules/nf-core/emmtyper/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'emmtyper'] )
}

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_emmtyper = Channel.empty()

    EMMTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(EMMTYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        EMMTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'emmtyper'], tsv]}.set{ ch_merge_emmtyper }
        CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')
        ch_merged_emmtyper = ch_merged_emmtyper.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = ch_merged_emmtyper
    versions = ch_versions
}
