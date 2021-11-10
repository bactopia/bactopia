//
// agrvate - Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
//
argvate_args = params.typing_only ? '--typing_only' : ''

include { AGRVATE as AGRVATE_MODULE } from '../../../modules/nf-core/modules/agrvate/main' addParams( options: [args: "${argvate_args}", is_module: true] )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow AGRVATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    AGRVATE_MODULE(fasta)
    ch_versions = ch_versions.mix(AGRVATE_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        AGRVATE_MODULE.out.summary.collect{meta, summary -> summary}.map{ summary -> [[id:'agrvate'], summary]}.set{ ch_merge_agrvate }
        CSVTK_CONCAT(ch_merge_agrvate, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }

    emit:
    tsv = AGRVATE_MODULE.out.summary
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
