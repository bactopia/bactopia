//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'staphopiasccmec')
options.args = params.hamming ? '--hamming' : ''

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE} from '../../../modules/nf-core/staphopiasccmec/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'staphopiasccmec'] )

workflow STAPHOPIASCCMEC {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    STAPHOPIASCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(STAPHOPIASCCMEC_MODULE.out.versions.first())

    // Merge results
    STAPHOPIASCCMEC_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'staphopiasccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
