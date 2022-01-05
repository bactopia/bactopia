//
// meningotype - Serotyping of Neisseria meningitidis
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'meningotype')
options.is_module = params.wf == 'meningotype' ? true : false
options.args = [
    params.finetype ? "--finetype" : "",
    params.porB ? "--porB" : "",
    params.bast ? "--bast" : "",
    params.mlst ? "--mlst" : "",
    params.all ? "--all" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../../modules/nf-core/modules/meningotype/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'meningotype'] )
}

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_meningotype = Channel.empty()

    MENINGOTYPE_MODULE(fasta)
    ch_versions = ch_versions.mix(MENINGOTYPE_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        MENINGOTYPE_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'meningotype'], tsv]}.set{ ch_merge_meningotype }
        CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')
        ch_merged_meningotype = ch_merged_meningotype.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = ch_merged_meningotype
    versions = ch_versions
}
