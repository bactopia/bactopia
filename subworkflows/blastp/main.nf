//
// blastp - Search against protein BLAST databases using protein queries
//
nextflow.preview.types = true

include { BLAST_BLASTP as BLASTP_MODULE } from '../../modules/blast/blastp/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTP {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    query // channel: [ fasta ]

    main:
    BLASTP_MODULE(fasta, query)

    // Merge results
    ch_merge_blastp = BLASTP_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'blastp'], tsv]}
    CSVTK_CONCAT(ch_merge_blastp, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BLASTP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = BLASTP_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = BLASTP_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = BLASTP_MODULE.out.nf_begin.mix(
        BLASTP_MODULE.out.nf_err,
        BLASTP_MODULE.out.nf_log,
        BLASTP_MODULE.out.nf_out,
        BLASTP_MODULE.out.nf_run,
        BLASTP_MODULE.out.nf_sh,
        BLASTP_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = BLASTP_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
