//
// gamma - A tool for identification, classification, and annotation of translated gene matches
//
include { GAMMA as GAMMA_MODULE } from '../../modules/gamma/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow GAMMA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    db

    main:
    GAMMA_MODULE(fasta, db)

    // Merge results
    GAMMA_MODULE.out.gamma.collect{_meta, gamma -> gamma}.map{ gamma -> [[id:'gamma'], gamma]}.set{ ch_merge_gamma }
    CSVTK_CONCAT(ch_merge_gamma, 'tsv', 'tsv')

    emit:
    // Individual outputs
    gamma = GAMMA_MODULE.out.gamma
    merged_gamma = CSVTK_CONCAT.out.csv
    psl = GAMMA_MODULE.out.psl
    fasta = GAMMA_MODULE.out.fasta
    gff = GAMMA_MODULE.out.gff

    // Generic aggregate outputs
    results = GAMMA_MODULE.out.gamma.mix(
        CSVTK_CONCAT.out.csv,
        GAMMA_MODULE.out.psl,
        GAMMA_MODULE.out.fasta,
        GAMMA_MODULE.out.gff
    )
    logs = GAMMA_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = GAMMA_MODULE.out.nf_begin.mix(
        GAMMA_MODULE.out.nf_err,
        GAMMA_MODULE.out.nf_log,
        GAMMA_MODULE.out.nf_out,
        GAMMA_MODULE.out.nf_run,
        GAMMA_MODULE.out.nf_sh,
        GAMMA_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = GAMMA_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
