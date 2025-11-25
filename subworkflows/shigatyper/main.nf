//
// shigatyper - Shigella serotype from Illumina or Oxford Nanopore reads
//
nextflow.preview.types = true

include { SHIGATYPER as SHIGATYPER_MODULE } from '../../modules/shigatyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SHIGATYPER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    SHIGATYPER_MODULE(reads)

    // Merge results
    ch_merge_shigatyper = SHIGATYPER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'shigatyper'], tsv]}
    CSVTK_CONCAT(ch_merge_shigatyper, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SHIGATYPER_MODULE.out.tsv
    hits = SHIGATYPER_MODULE.out.hits
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = SHIGATYPER_MODULE.out.tsv.mix(
        SHIGATYPER_MODULE.out.hits,
        CSVTK_CONCAT.out.csv
    )
    logs = SHIGATYPER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SHIGATYPER_MODULE.out.nf_begin.mix(
        SHIGATYPER_MODULE.out.nf_err,
        SHIGATYPER_MODULE.out.nf_log,
        SHIGATYPER_MODULE.out.nf_out,
        SHIGATYPER_MODULE.out.nf_run,
        SHIGATYPER_MODULE.out.nf_sh,
        SHIGATYPER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SHIGATYPER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
