//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'staphopiasccmec')
options.is_module = params.wf == 'staphopiasccmec' ? true : false
options.args = params.hamming ? '--hamming' : ''

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE} from '../../../modules/nf-core/modules/staphopiasccmec/main' addParams( options: options )
if (params.is_subworkflow) include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'staphopiasccmec'] )

workflow STAPHOPIASCCMEC {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_staphopiasccmec = Channel.empty()

    STAPHOPIASCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(STAPHOPIASCCMEC_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        STAPHOPIASCCMEC_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'staphopiasccmec'], tsv]}.set{ ch_merge_sccmec }
        CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
        ch_merged_staphopiasccmec = ch_merged_staphopiasccmec.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv = ch_merged_staphopiasccmec
    versions = ch_versions
}
