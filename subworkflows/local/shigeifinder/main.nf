//
// shigeifinder - Shigella and EIEC serotyping from assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigeifinder')

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../../modules/nf-core/shigeifinder/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'shigeifinder'] )

workflow SHIGEIFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SHIGEIFINDER_MODULE(fasta)
    ch_versions = ch_versions.mix(SHIGEIFINDER_MODULE.out.versions.first())

    // Merge results
    SHIGEIFINDER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'shigeifinder'], tsv]}.set{ ch_merge_shigeifinder }
    CSVTK_CONCAT(ch_merge_shigeifinder, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SHIGEIFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
