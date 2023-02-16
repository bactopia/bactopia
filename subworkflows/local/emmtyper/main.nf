//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'emmtyper')
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
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'emmtyper'] )

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    EMMTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(EMMTYPER_MODULE.out.versions.first())

    // Merge results
    EMMTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'emmtyper'], tsv]}.set{ ch_merge_emmtyper }
    CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
