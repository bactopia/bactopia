//
// blastn - Search against nucleotide BLAST databases using nucleotide queries
//
include { BLAST_BLASTN as BLASTN_MODULE } from '../../modules/blast/blastn/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTN {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    query // channel: [ fasta ]

    main:
    BLASTN_MODULE(fasta, query)

    // Merge results
    BLASTN_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'blastn'], tsv]}.set{ ch_merge_blastn }
    CSVTK_CONCAT(ch_merge_blastn, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BLASTN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = BLASTN_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = BLASTN_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = BLASTN_MODULE.out.nf_begin.mix(
        BLASTN_MODULE.out.nf_err,
        BLASTN_MODULE.out.nf_log,
        BLASTN_MODULE.out.nf_out,
        BLASTN_MODULE.out.nf_run,
        BLASTN_MODULE.out.nf_sh,
        BLASTN_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = BLASTN_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
