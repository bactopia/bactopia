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
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    DEFENSEFINDER_UPDATE()
    ch_versions = ch_versions.mix(DEFENSEFINDER_UPDATE.out.versions)
    ch_logs = ch_logs.mix(DEFENSEFINDER_UPDATE.out.logs)
    ch_versions = ch_versions.mix(DEFENSEFINDER_RUN.out.versions)
    ch_logs = ch_logs.mix(DEFENSEFINDER_RUN.out.logs)

    ch_versions = ch_versions.mix(GENES_CONCAT.out.versions)
    ch_logs = ch_logs.mix(GENES_CONCAT.out.logs)

    ch_versions = ch_versions.mix(HMMER_CONCAT.out.versions)
    ch_logs = ch_logs.mix(HMMER_CONCAT.out.logs)

    ch_versions = ch_versions.mix(SYSTEMS_CONCAT.out.versions)
    ch_logs = ch_logs.mix(SYSTEMS_CONCAT.out.logs)

    emit:
    ch_merged_genes = ch_merged_genes.mix(GENES_CONCAT.out.csv)
    ch_merged_hmmer = ch_merged_hmmer.mix(HMMER_CONCAT.out.csv)
    ch_merged_systems = ch_merged_systems.mix(SYSTEMS_CONCAT.out.csv)
    genes_tsv = DEFENSEFINDER_RUN.out.genes_tsv
    hmmer_tsv = DEFENSEFINDER_RUN.out.genes_tsv
    macsydata_raw = DEFENSEFINDER_RUN.out.macsydata_raw
    merged_genes_tsv = ch_merged_genes
    merged_hmmer_tsv = ch_merged_hmmer
    merged_systems_tsv = ch_merged_systems
    proteins = DEFENSEFINDER_RUN.out.proteins
    proteins_index = DEFENSEFINDER_RUN.out.proteins_index
    systems_tsv = DEFENSEFINDER_RUN.out.systems_tsv
    logs = ch_logs
    nf_logs = DEFENSEFINDER_UPDATE.out.nf_begin.mix(
        DEFENSEFINDER_UPDATE.out.nf_err,
        DEFENSEFINDER_UPDATE.out.nf_log,
        DEFENSEFINDER_UPDATE.out.nf_out,
        DEFENSEFINDER_UPDATE.out.nf_run,
        DEFENSEFINDER_UPDATE.out.nf_sh,
        DEFENSEFINDER_UPDATE.out.nf_trace
    )
    versions = ch_versions
}
