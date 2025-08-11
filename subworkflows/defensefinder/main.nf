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
    ch_nf_logs = Channel.empty()
    ch_merged_genes = Channel.empty()
    ch_merged_hmmer = Channel.empty()
    ch_merged_systems = Channel.empty()

    DEFENSEFINDER_UPDATE()
    ch_versions = ch_versions.mix(DEFENSEFINDER_UPDATE.out.versions)
    ch_logs = ch_logs.mix(DEFENSEFINDER_UPDATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(DEFENSEFINDER_UPDATE.out.nf_logs)

    DEFENSEFINDER_RUN(fasta, DEFENSEFINDER_UPDATE.out.db)
    ch_versions = ch_versions.mix(DEFENSEFINDER_RUN.out.versions)
    ch_logs = ch_logs.mix(DEFENSEFINDER_RUN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(DEFENSEFINDER_RUN.out.nf_logs)

    // Merge results
    DEFENSEFINDER_RUN.out.genes_tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'defensefinder-genes'], tsv] }.set{ ch_merge_genes }
    GENES_CONCAT(ch_merge_genes, 'tsv', 'tsv')
    ch_merged_genes = ch_merged_genes.mix(GENES_CONCAT.out.csv)
    ch_versions = ch_versions.mix(GENES_CONCAT.out.versions)
    ch_logs = ch_logs.mix(GENES_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GENES_CONCAT.out.nf_logs)

    DEFENSEFINDER_RUN.out.hmmer_tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'defensefinder-hmmer'], tsv] }.set{ ch_merge_hmmer }
    HMMER_CONCAT(ch_merge_hmmer, 'tsv', 'tsv')
    ch_merged_hmmer = ch_merged_hmmer.mix(HMMER_CONCAT.out.csv)
    ch_versions = ch_versions.mix(HMMER_CONCAT.out.versions)
    ch_logs = ch_logs.mix(HMMER_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(HMMER_CONCAT.out.nf_logs)

    DEFENSEFINDER_RUN.out.systems_tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'defensefinder-systems'], tsv] }.set{ ch_merge_systems }
    SYSTEMS_CONCAT(ch_merge_systems, 'tsv', 'tsv')
    ch_merged_systems = ch_merged_systems.mix(SYSTEMS_CONCAT.out.csv)
    ch_versions = ch_versions.mix(SYSTEMS_CONCAT.out.versions)
    ch_logs = ch_logs.mix(SYSTEMS_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SYSTEMS_CONCAT.out.nf_logs)

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
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
