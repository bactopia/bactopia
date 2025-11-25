//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
nextflow.preview.types = true

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow HPSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    HPSUISSERO_MODULE(fasta)

    // Merge results
    ch_merge_hpsuissero = HPSUISSERO_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'hpsuissero'], tsv]}
    CSVTK_CONCAT(ch_merge_hpsuissero, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = HPSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = HPSUISSERO_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = HPSUISSERO_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        HPSUISSERO_MODULE.out.nf_begin,
        HPSUISSERO_MODULE.out.nf_err,
        HPSUISSERO_MODULE.out.nf_log,
        HPSUISSERO_MODULE.out.nf_out,
        HPSUISSERO_MODULE.out.nf_run,
        HPSUISSERO_MODULE.out.nf_sh,
        HPSUISSERO_MODULE.out.nf_trace
    )
    versions = HPSUISSERO_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
