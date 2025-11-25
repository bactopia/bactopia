//
// agrvate - Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
//
nextflow.preview.types = true

include { AGRVATE as AGRVATE_MODULE } from '../../modules/agrvate/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow AGRVATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    AGRVATE_MODULE(fasta)

    // Merge results
    ch_merge_agrvate = AGRVATE_MODULE.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:'agrvate'], summary]}
    CSVTK_CONCAT(ch_merge_agrvate, 'tsv', 'tsv')

    emit:
    // Individual output
    tsv = AGRVATE_MODULE.out.summary
    supplemental = AGRVATE_MODULE.out.supplemental
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate output
    results = AGRVATE_MODULE.out.summary.mix(
        AGRVATE_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    )
    logs = AGRVATE_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = AGRVATE_MODULE.out.nf_begin.mix(
        AGRVATE_MODULE.out.nf_err,
        AGRVATE_MODULE.out.nf_log,
        AGRVATE_MODULE.out.nf_out,
        AGRVATE_MODULE.out.nf_run,
        AGRVATE_MODULE.out.nf_sh,
        AGRVATE_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = AGRVATE_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
