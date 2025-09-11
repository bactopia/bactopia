//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE} from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow STAPHOPIASCCMEC {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    STAPHOPIASCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(STAPHOPIASCCMEC_MODULE.out.versions)

    // Merge results
    STAPHOPIASCCMEC_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'staphopiasccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = STAPHOPIASCCMEC_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = STAPHOPIASCCMEC_MODULE.out.nf_begin.mix(
        STAPHOPIASCCMEC_MODULE.out.nf_err,
        STAPHOPIASCCMEC_MODULE.out.nf_log,
        STAPHOPIASCCMEC_MODULE.out.nf_out,
        STAPHOPIASCCMEC_MODULE.out.nf_run,
        STAPHOPIASCCMEC_MODULE.out.nf_sh,
        STAPHOPIASCCMEC_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions
}
