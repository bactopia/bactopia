//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
emmtyper_opts = [
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

include { EMMTYPER as EMMTYPER_MODULE } from '../../../modules/nf-core/modules/emmtyper/main' addParams( options: [ args: "${emmtyper_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    EMMTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(EMMTYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        EMMTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'emmtyper'], tsv]}.set{ ch_merge_emmtyper }
        CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')
    }

    emit:
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
