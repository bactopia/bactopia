//
// seqsero2 - Salmonella serotype prediction from reads or assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'seqsero2')
options.is_module = params.wf == 'seqsero2' ? true : false
options.args = [
    "-m ${params.run_mode}",
    params.input_type == "assembly" ? "-t 4" : "-t 2",
    "-b ${params.bwa_mode}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
include { SEQSERO2 as SEQSERO2_MODULE } from '../../../modules/nf-core/modules/seqsero2/main' addParams( options: options )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'seqsero2'] )
}

workflow SEQSERO2 {
    take:
    seqs // channel: [ val(meta), [ fastqs or assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_seqsero2 = Channel.empty()

    SEQSERO2_MODULE(seqs)
    ch_versions = ch_versions.mix(SEQSERO2_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SEQSERO2_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'seqsero2'], tsv]}.set{ ch_merge_seqsero2 }
        CSVTK_CONCAT(ch_merge_seqsero2, 'tsv', 'tsv')
        ch_merged_seqsero2 = ch_merged_seqsero2.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SEQSERO2_MODULE.out.tsv
    txt = SEQSERO2_MODULE.out.txt
    merged_tsv = ch_merged_seqsero2
    versions = ch_versions
}
