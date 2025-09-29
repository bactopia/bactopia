//
// ngmaster - Multi-antigen sequence typing for Neisseria gonorrhoeae
//
include { NGMASTER as NGMASTER_MODULE } from '../../modules/ngmaster/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow NGMASTER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    NGMASTER_MODULE(fasta)

    // Merge results
    NGMASTER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'ngmaster'], tsv]}.set{ ch_merge_ngmaster }
    CSVTK_CONCAT(ch_merge_ngmaster, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = NGMASTER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = NGMASTER_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = NGMASTER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = NGMASTER_MODULE.out.nf_begin.mix(
        NGMASTER_MODULE.out.nf_err,
        NGMASTER_MODULE.out.nf_log,
        NGMASTER_MODULE.out.nf_out,
        NGMASTER_MODULE.out.nf_run,
        NGMASTER_MODULE.out.nf_sh,
        NGMASTER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = NGMASTER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
