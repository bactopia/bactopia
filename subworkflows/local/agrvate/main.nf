//
// agrvate - Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'agrvate')
options.args = params.typing_only ? '--typing_only' : ''
options.subdir = "${params.run_name}"

include { AGRVATE as AGRVATE_MODULE } from '../../../modules/nf-core/agrvate/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [args: '-C "$"', logs_subdir: 'agrvate-concat', process_name: params.merge_folder] )

workflow AGRVATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_agrvate = Channel.empty()

    AGRVATE_MODULE(fasta)
    ch_versions = ch_versions.mix(AGRVATE_MODULE.out.versions.first())

    AGRVATE_MODULE.out.summary.collect{meta, summary -> summary}.map{ summary -> [[id:'agrvate'], summary]}.set{ ch_merge_agrvate }
    CSVTK_CONCAT(ch_merge_agrvate, 'tsv', 'tsv')
    ch_merged_agrvate = ch_merged_agrvate.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = AGRVATE_MODULE.out.summary
    merged_tsv = ch_merged_agrvate
    versions = ch_versions
}
