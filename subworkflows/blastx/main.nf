//
// blastx - Search against protein BLAST databases using translated nucleotide queries
//
nextflow.preview.types = true

include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTX {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    query // channel: [ fasta ]

    main:
    BLASTX_MODULE(fasta, query)

    // Merge results
    ch_merge_blastx = BLASTX_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'blastx'], tsv]}
    CSVTK_CONCAT(ch_merge_blastx, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = BLASTX_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = BLASTX_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = BLASTX_MODULE.out.nf_begin.mix(
        BLASTX_MODULE.out.nf_err,
        BLASTX_MODULE.out.nf_log,
        BLASTX_MODULE.out.nf_out,
        BLASTX_MODULE.out.nf_run,
        BLASTX_MODULE.out.nf_sh,
        BLASTX_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = BLASTX_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
