//
// sistr - Serovar prediction of Salmonella assemblies
//
include { SISTR as SISTR_MODULE } from '../../modules/sistr/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SISTR {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    SISTR_MODULE(fasta)

    // Merge results
    SISTR_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'sistr'], tsv]}.set{ ch_merge_sistr }
    CSVTK_CONCAT(ch_merge_sistr, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SISTR_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    allele_fasta = SISTR_MODULE.out.allele_fasta
    allele_json = SISTR_MODULE.out.allele_json
    cgmlst_csv = SISTR_MODULE.out.cgmlst_csv

    // Generic aggregate outputs
    results = SISTR_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv,
        SISTR_MODULE.out.allele_fasta,
        SISTR_MODULE.out.allele_json,
        SISTR_MODULE.out.cgmlst_csv
    )
    logs = SISTR_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SISTR_MODULE.out.nf_begin.mix(
        SISTR_MODULE.out.nf_err,
        SISTR_MODULE.out.nf_log,
        SISTR_MODULE.out.nf_out,
        SISTR_MODULE.out.nf_run,
        SISTR_MODULE.out.nf_sh,
        SISTR_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SISTR_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
