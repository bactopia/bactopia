//
// gamma - A tool for identification, classification, and annotation of translated gene matches
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'gamma')
options.is_module = params.wf == 'gamma' ? true : false
options.args = [
    params.write_fasta ? "--fasta" : "",
    params.write_fasta ? "--gff" : "",
    params.extended ? "--extended" : "",
    params.all_matches ? "--all" : "",
    "--percent_identity ${params.percent_identity}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
GAMMA_DB = params.gamma_db ? file(params.gamma_db) : []

include { GAMMA as GAMMA_MODULE } from '../../../modules/nf-core/modules/gamma/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'gamma'] )
}

workflow GAMMA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_gamma = Channel.empty()

    GAMMA_MODULE(fasta, GAMMA_DB)
    ch_versions = ch_versions.mix(GAMMA_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        GAMMA_MODULE.out.gamma.collect{meta, gamma -> gamma}.map{ gamma -> [[id:'gamma'], gamma]}.set{ ch_merge_gamma }
        CSVTK_CONCAT(ch_merge_gamma, 'tsv', 'tsv')
        ch_merged_gamma = ch_merged_gamma.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }


    emit:
    gamma = GAMMA_MODULE.out.gamma
    merged_gamma = ch_merged_gamma
    psl = GAMMA_MODULE.out.psl
    fasta = GAMMA_MODULE.out.fasta
    gff = GAMMA_MODULE.out.gff
    versions = ch_versions // channel: [ versions.yml ]
}
