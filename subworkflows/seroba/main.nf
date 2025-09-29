//
// seroba - Serotyping of Streptococcus pneumoniae from sequence reads
//
include { SEROBA_RUN } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SEROBA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    SEROBA_RUN(fasta)

    // Merge results
    SEROBA_RUN.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'seroba'], tsv]}.set{ ch_merge_seroba }
    CSVTK_CONCAT(ch_merge_seroba, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SEROBA_RUN.out.tsv
    txt = SEROBA_RUN.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = SEROBA_RUN.out.tsv.mix(
        SEROBA_RUN.out.txt,
        CSVTK_CONCAT.out.csv
    )
    logs = SEROBA_RUN.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SEROBA_RUN.out.nf_begin.mix(
        SEROBA_RUN.out.nf_err,
        SEROBA_RUN.out.nf_log,
        SEROBA_RUN.out.nf_out,
        SEROBA_RUN.out.nf_run,
        SEROBA_RUN.out.nf_sh,
        SEROBA_RUN.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SEROBA_RUN.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
