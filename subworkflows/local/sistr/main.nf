//
// sistr - Serovar prediction of Salmonella assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'sistr')
options.is_module = params.wf == 'sistr' ? true : false
options.args = [
    params.full_cgmlst ? "--use-full-cgmlst-db" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { SISTR as SISTR_MODULE } from '../../../modules/nf-core/modules/sistr/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'sistr'] )
}

workflow SISTR {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_sistr = Channel.empty()

    SISTR_MODULE(fasta)
    ch_versions = ch_versions.mix(SISTR_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SISTR_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'sistr'], tsv]}.set{ ch_merge_sistr }
        CSVTK_CONCAT(ch_merge_sistr, 'tsv', 'tsv')
        ch_merged_sistr = ch_merged_sistr.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SISTR_MODULE.out.tsv
    merged_tsv = ch_merged_sistr
    allele_fasta = SISTR_MODULE.out.allele_fasta
    allele_json = SISTR_MODULE.out.allele_json
    cgmlst_csv = SISTR_MODULE.out.cgmlst_csv
    versions = ch_versions
}
