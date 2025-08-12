//
// gtdb - Identify marker genes and assign taxonomic classifications
//
include { GTDBTK_SETUPDB as SETUPDB } from '../../modules/gtdbtk/setupdb/main'
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow GTDB {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    DATABASE = ! params.download_gtdb ? file(params.gtdb) : []

    if (params.download_gtdb) {
        // Force CLASSIFY to wait
        SETUPDB()
        ch_versions = ch_versions.mix(SETUPDB.out.versions)
        ch_logs = ch_logs.mix(SETUPDB.out.logs)

        if (params.gtdb_save_as_tarball) {
            CLASSIFY(fasta, SETUPDB.out.db_tarball)
        } else {
            CLASSIFY(fasta, SETUPDB.out.db)
        }
    } else {
        CLASSIFY(fasta, DATABASE)
    }
    
    ch_versions = ch_versions.mix(CLASSIFY.out.versions.first())
    ch_logs = ch_logs.mix(CLASSIFY.out.logs)
    
    // Merge results
    CLASSIFY.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'gtdb'], tsv]}.set{ ch_merge_gtdb }
    CSVTK_CONCAT(ch_merge_gtdb, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = CLASSIFY.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    results = CLASSIFY.out.results
    logs = ch_logs
    nf_logs = CLASSIFY.out.nf_begin.mix(
        CLASSIFY.out.nf_err,
        CLASSIFY.out.nf_log,
        CLASSIFY.out.nf_out,
        CLASSIFY.out.nf_run,
        CLASSIFY.out.nf_sh,
        CLASSIFY.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions
}
