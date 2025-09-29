//
// fastani - fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
//
include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow FASTANI {
    take:
    query // channel: [ val(meta), [ fasta ] ]
    reference // channel: [ val(meta), [ fasta ] ]

    main:
    query.collect{_meta, fasta -> fasta}.map{ fasta -> [[id:'query'], fasta]}.set{ ch_fastani_query }
    reference.map{_meta, fasta -> fasta}.set{ ch_fastani_reference }
    FASTANI_MODULE(ch_fastani_query, ch_fastani_reference)

    // Merge results
    FASTANI_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'fastani'], tsv]}.set{ ch_merge_fastani }
    CSVTK_CONCAT(ch_merge_fastani, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = FASTANI_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = FASTANI_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = FASTANI_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = FASTANI_MODULE.out.nf_begin.mix(
        FASTANI_MODULE.out.nf_err,
        FASTANI_MODULE.out.nf_log,
        FASTANI_MODULE.out.nf_out,
        FASTANI_MODULE.out.nf_run,
        FASTANI_MODULE.out.nf_sh,
        FASTANI_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = FASTANI_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
