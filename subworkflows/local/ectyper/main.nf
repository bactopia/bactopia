//
// ectyper - In-silico prediction of Escherichia coli serotype
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ectyper')
options.is_module = params.wf == 'ectyper' ? true : false
options.args = [
    params.verify ? "--verify" : "",
    params.print_alleles  ? "-s" : "",
    "-opid ${params.opid}",
    "-opcov ${params.opcov}",
    "-hpid ${params.hpid}",
    "-hpcov ${params.hpcov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ECTYPER as ECTYPER_MODULE } from '../../../modules/nf-core/modules/ectyper/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'ectyper'] )
}

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_ectyper = Channel.empty()

    ECTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(ECTYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        ECTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ectyper'], tsv]}.set{ ch_merge_ectyper }
        CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')
        ch_merged_ectyper = ch_merged_ectyper.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = ECTYPER_MODULE.out.tsv
    merged_tsv = ch_merged_ectyper
    txt = ECTYPER_MODULE.out.txt
    versions = ch_versions
}
