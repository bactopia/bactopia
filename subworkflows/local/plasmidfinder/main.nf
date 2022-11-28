//
// plasmidfinder - Plasmid identification from assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'plasmidfinder')
options.is_module = params.wf == 'plasmidfinder' ? true : false
options.args = [
    "-l ${params.pf_mincov}",
    "-t ${params.pf_threshold}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../../modules/nf-core/plasmidfinder/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'plasmidfinder'] )
}

workflow PLASMIDFINDER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_plasmidfinder = Channel.empty()

    PLASMIDFINDER_MODULE(fasta)
    ch_versions = ch_versions.mix(PLASMIDFINDER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        PLASMIDFINDER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'plasmidfinder'], tsv]}.set{ ch_merge_plasmidfinder }
        CSVTK_CONCAT(ch_merge_plasmidfinder, 'tsv', 'tsv')
        ch_merged_plasmidfinder = ch_merged_plasmidfinder.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    json = PLASMIDFINDER_MODULE.out.json
    txt = PLASMIDFINDER_MODULE.out.txt
    tsv = PLASMIDFINDER_MODULE.out.tsv
    merged_tsv = ch_merged_plasmidfinder
    genome_seq = PLASMIDFINDER_MODULE.out.genome_seq
    plasmid_seq = PLASMIDFINDER_MODULE.out.plasmid_seq
    versions = ch_versions
}
