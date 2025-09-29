//
// gtdb - Identify marker genes and assign taxonomic classifications
//
include { GTDBTK_DOWNLOAD as DOWNLOAD } from '../../modules/gtdbtk/download/main'
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow GTDB {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]
    database
    download_gtdb
    save_as_tarball

    main:
    if (download_gtdb) {
        // Force CLASSIFY to wait
        DOWNLOAD()

        if (save_as_tarball) {
            CLASSIFY(fasta, DOWNLOAD.out.db_tarball)
        } else {
            CLASSIFY(fasta, DOWNLOAD.out.db)
        }
    } else {
        CLASSIFY(fasta, database)
    }
    
    // Merge results
    CLASSIFY.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'gtdb'], tsv]}.set{ ch_merge_gtdb }
    CSVTK_CONCAT(ch_merge_gtdb, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = CLASSIFY.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = CLASSIFY.out.tsv.mix(
        CLASSIFY.out.supplemental,
        CSVTK_CONCAT.out.csv
    )
    logs = CLASSIFY.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = CLASSIFY.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
