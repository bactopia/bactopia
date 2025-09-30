//
// defensefinder - Systematic search of all known anti-phage systems
//
include { DEFENSEFINDER_UPDATE } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN    } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'

workflow DEFENSEFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(fasta, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    DEFENSEFINDER_RUN.out.genes_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-genes'], tsv]}.set{ ch_merge_genes }
    GENES_CONCAT(ch_merge_genes, 'tsv', 'tsv')

    DEFENSEFINDER_RUN.out.hmmer_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-hmmer'], tsv]}.set{ ch_merge_hmmer }
    HMMER_CONCAT(ch_merge_hmmer, 'tsv', 'tsv')

    DEFENSEFINDER_RUN.out.systems_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-systems'], tsv]}.set{ ch_merge_systems }
    SYSTEMS_CONCAT(ch_merge_systems, 'tsv', 'tsv')

    emit:
    // Individual outputs
    genes_tsv = DEFENSEFINDER_RUN.out.genes_tsv
    merged_genes_tsv = GENES_CONCAT.out.csv
    hmmer_tsv = DEFENSEFINDER_RUN.out.hmmer_tsv
    merged_hmmer_tsv = HMMER_CONCAT.out.csv
    systems_tsv = DEFENSEFINDER_RUN.out.systems_tsv
    merged_systems_tsv = SYSTEMS_CONCAT.out.csv
    proteins = DEFENSEFINDER_RUN.out.proteins
    proteins_index = DEFENSEFINDER_RUN.out.proteins_index
    macsydata_raw = DEFENSEFINDER_RUN.out.macsydata_raw

    // Generic aggregate outputs
    results = DEFENSEFINDER_RUN.out.genes_tsv.mix(
        GENES_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.hmmer_tsv,
        HMMER_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.systems_tsv,
        SYSTEMS_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.proteins,
        DEFENSEFINDER_RUN.out.proteins_index,
        DEFENSEFINDER_RUN.out.macsydata_raw
    )
    logs = DEFENSEFINDER_RUN.out.logs.mix(
        GENES_CONCAT.out.logs,
        HMMER_CONCAT.out.logs,
        SYSTEMS_CONCAT.out.logs
    )
    nf_logs = DEFENSEFINDER_RUN.out.nf_begin.mix(
        DEFENSEFINDER_RUN.out.nf_err,
        DEFENSEFINDER_RUN.out.nf_log,
        DEFENSEFINDER_RUN.out.nf_out,
        DEFENSEFINDER_RUN.out.nf_run,
        DEFENSEFINDER_RUN.out.nf_sh,
        DEFENSEFINDER_RUN.out.nf_trace,
        GENES_CONCAT.out.nf_begin,
        GENES_CONCAT.out.nf_err,
        GENES_CONCAT.out.nf_log,
        GENES_CONCAT.out.nf_out,
        GENES_CONCAT.out.nf_run,
        GENES_CONCAT.out.nf_sh,
        GENES_CONCAT.out.nf_trace,
        HMMER_CONCAT.out.nf_begin,
        HMMER_CONCAT.out.nf_err,
        HMMER_CONCAT.out.nf_log,
        HMMER_CONCAT.out.nf_out,
        HMMER_CONCAT.out.nf_run,
        HMMER_CONCAT.out.nf_sh,
        HMMER_CONCAT.out.nf_trace,
        SYSTEMS_CONCAT.out.nf_begin,
        SYSTEMS_CONCAT.out.nf_err,
        SYSTEMS_CONCAT.out.nf_log,
        SYSTEMS_CONCAT.out.nf_out,
        SYSTEMS_CONCAT.out.nf_run,
        SYSTEMS_CONCAT.out.nf_sh,
        SYSTEMS_CONCAT.out.nf_trace
    )
    versions = DEFENSEFINDER_RUN.out.versions.mix(
        GENES_CONCAT.out.versions,
        HMMER_CONCAT.out.versions,
        SYSTEMS_CONCAT.out.versions
    )
}
