//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
nextflow.preview.types = true

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE} from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow STAPHOPIASCCMEC {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    STAPHOPIASCCMEC_MODULE(fasta)

    // Merge results
    ch_merge_sccmec = STAPHOPIASCCMEC_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'staphopiasccmec'], tsv]}
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = STAPHOPIASCCMEC_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
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
    versions = STAPHOPIASCCMEC_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
