//
// defensefinder - Systematic search of all known anti-phage systems
//
include { DEFENSEFINDER_UPDATE } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'

workflow DEFENSEFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_logs = Channel.empty()
    ch_versions = Channel.empty()
    ch_merged_genes = Channel.empty()
    ch_merged_hmmer = Channel.empty()
    ch_merged_systems = Channel.empty()

    // Update Defensefinder
    DEFENSEFINDER_UPDATE()

    // Run Defensefinder
    DEFENSEFINDER_RUN(fasta, DEFENSEFINDER_UPDATE.out.db)
    ch_versions = ch_versions.mix(DEFENSEFINDER_RUN.out.versions)
    ch_logs = ch_logs.mix(DEFENSEFINDER_RUN.out.logs)

    // Merge results
    DEFENSEFINDER_RUN.out.genes_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-genes'], tsv]}.set{ ch_merge_genes }
    GENES_CONCAT(ch_merge_genes, 'tsv', 'tsv')
    ch_merged_genes = ch_merged_genes.mix(GENES_CONCAT.out.csv)
    ch_logs = ch_logs.mix(GENES_CONCAT.out.logs)
    ch_versions = ch_versions.mix(GENES_CONCAT.out.versions)

    DEFENSEFINDER_RUN.out.hmmer_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-hmmer'], tsv]}.set{ ch_merge_hmmer }
    HMMER_CONCAT(ch_merge_hmmer, 'tsv', 'tsv')
    ch_merged_hmmer = ch_merged_hmmer.mix(HMMER_CONCAT.out.csv)
    ch_logs = ch_logs.mix(HMMER_CONCAT.out.logs)
    ch_versions = ch_versions.mix(HMMER_CONCAT.out.versions)

    DEFENSEFINDER_RUN.out.systems_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-systems'], tsv]}.set{ ch_merge_systems }
    SYSTEMS_CONCAT(ch_merge_systems, 'tsv', 'tsv')
    ch_merged_systems = ch_merged_systems.mix(SYSTEMS_CONCAT.out.csv)
    ch_logs = ch_logs.mix(SYSTEMS_CONCAT.out.logs)
    ch_versions = ch_versions.mix(SYSTEMS_CONCAT.out.versions)

    emit:
    genes_tsv = DEFENSEFINDER_RUN.out.genes_tsv
    merged_genes_tsv = ch_merged_genes
    hmmer_tsv = DEFENSEFINDER_RUN.out.genes_tsv
    merged_hmmer_tsv = ch_merged_hmmer
    systems_tsv = DEFENSEFINDER_RUN.out.systems_tsv
    merged_systems_tsv = ch_merged_systems
    proteins = DEFENSEFINDER_RUN.out.proteins
    proteins_index = DEFENSEFINDER_RUN.out.proteins_index
    macsydata_raw = DEFENSEFINDER_RUN.out.macsydata_raw
    logs = ch_logs
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
    versions = ch_versions
}
