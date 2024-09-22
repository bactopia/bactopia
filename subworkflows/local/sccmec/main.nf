//
// sccmec - A tool for typing SCCmec cassettes in assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'sccmec')
options.args = [
    "--min-targets-pident ${params.sccmec_min_targets_pident}",
    "--min-targets-coverage ${params.sccmec_min_targets_coverage}",
    "--min-regions-pident ${params.sccmec_min_regions_pident}",
    "--min-regions-coverage ${params.sccmec_min_regions_coverage}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { SCCMEC as SCCMEC_MODULE } from '../../../modules/nf-core/sccmec/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'sccmec-concat', process_name: params.merge_folder] )

workflow SCCMEC {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(SCCMEC_MODULE.out.versions.first())

    // Merge results
    SCCMEC_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'sccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    targets = SCCMEC_MODULE.out.targets
    target_details = SCCMEC_MODULE.out.target_details
    regions = SCCMEC_MODULE.out.regions
    regions_details = SCCMEC_MODULE.out.regions_details
    versions = ch_versions
}
