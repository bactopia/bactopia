//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
staphopiasccmec_args = params.hamming ? '--hamming' : ''

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE} from '../../../modules/nf-core/modules/staphopiasccmec/main' addParams( options: [args: "${staphopiasccmec_args}", is_module: true] )
if (params.is_subworkflow) include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )

workflow STAPHOPIASCCMEC {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    STAPHOPIASCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(STAPHOPIASCCMEC_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        STAPHOPIASCCMEC_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'staphopiasccmec'], tsv]}.set{ ch_merge_sccmec }
        CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }

    emit:
    tsv = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
