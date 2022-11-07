//
// pbptyper - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pbptyper')
options.is_module = params.wf == 'pbptyper' ? true : false
options.args = [
    "--min_pident ${params.min_pident}",
    "--min_coverage ${params.min_coverage}",
    "--min_ani ${params.min_ani}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PBPTYPER as PBPTYPER_MODULE } from '../../../modules/nf-core/pbptyper/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'pbptyper'] )
}

workflow PBPTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_pbptyper = Channel.empty()

    PBPTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(PBPTYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        PBPTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'pbptyper'], tsv]}.set{ ch_merge_pbptyper }
        CSVTK_CONCAT(ch_merge_pbptyper, 'tsv', 'tsv')
        ch_merged_pbptyper = ch_merged_pbptyper.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = PBPTYPER_MODULE.out.tsv
    merged_tsv = ch_merged_pbptyper
    blast = PBPTYPER_MODULE.out.blast
    versions = ch_versions
}
