//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
nextflow.preview.types = true

include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    blastdb

    main:
    EMMTYPER_MODULE(fasta, blastdb)

    // Merge results
    ch_merge_emmtyper = EMMTYPER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'emmtyper'], tsv]}
    CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = EMMTYPER_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = EMMTYPER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        EMMTYPER_MODULE.out.nf_begin,
        EMMTYPER_MODULE.out.nf_err,
        EMMTYPER_MODULE.out.nf_log,
        EMMTYPER_MODULE.out.nf_out,
        EMMTYPER_MODULE.out.nf_run,
        EMMTYPER_MODULE.out.nf_sh,
        EMMTYPER_MODULE.out.nf_trace
    )
    versions = EMMTYPER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
