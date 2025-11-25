//
// busco - Assembly completeness based on evolutionarily informed expectations
// 
include { BUSCO as BUSCO_MODULE } from '../../modules/busco/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BUSCO {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]
    busco_lineage

    main:
    BUSCO_MODULE(fasta)

    // Merge the results
    ch_merge_busco = BUSCO_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:"busco-${busco_lineage}"], tsv]}
    CSVTK_CONCAT(ch_merge_busco, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BUSCO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = BUSCO_MODULE.out.tsv.mix(
        BUSCO_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    )
    logs = BUSCO_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = BUSCO_MODULE.out.nf_begin.mix(
        BUSCO_MODULE.out.nf_err,
        BUSCO_MODULE.out.nf_log,
        BUSCO_MODULE.out.nf_out,
        BUSCO_MODULE.out.nf_run,
        BUSCO_MODULE.out.nf_sh,
        BUSCO_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = BUSCO_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
