//
// defensefinder - Systematic search of all known anti-phage systems
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'defensefinder')
options.args = [
    params.df_preserveraw ? "--preserve-raw" : "",
    params.df_nocutga ? "--no-cut-ga" : "",
    "--coverage ${params.df_coverage}",
    "--db-type ${params.df_dbtype}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { DEFENSEFINDER_UPDATE } from '../../../modules/nf-core/defensefinder/update/main' addParams()
include { DEFENSEFINDER_RUN } from '../../../modules/nf-core/defensefinder/run/main' addParams( options: options )
include { CSVTK_CONCAT as GENES_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'defensefinder-genes-concat', process_name: params.merge_folder] )
include { CSVTK_CONCAT as HMMER_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'defensefinder-hmmer-concat', process_name: params.merge_folder] )
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'defensefinder-systems-concat', process_name: params.merge_folder] )

workflow DEFENSEFINDER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_genes = Channel.empty()
    ch_merged_hmmer = Channel.empty()
    ch_merged_systems = Channel.empty()

    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(fasta, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    DEFENSEFINDER_RUN.out.genes_tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-genes'], tsv]}.set{ ch_merge_genes }
    GENES_CONCAT(ch_merge_genes, 'tsv', 'tsv')
    ch_merged_genes = ch_merged_genes.mix(GENES_CONCAT.out.csv)
    ch_versions = ch_versions.mix(GENES_CONCAT.out.versions)

    DEFENSEFINDER_RUN.out.hmmer_tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-hmmer'], tsv]}.set{ ch_merge_hmmer }
    HMMER_CONCAT(ch_merge_hmmer, 'tsv', 'tsv')
    ch_merged_hmmer = ch_merged_hmmer.mix(HMMER_CONCAT.out.csv)
    ch_versions = ch_versions.mix(HMMER_CONCAT.out.versions)

    DEFENSEFINDER_RUN.out.systems_tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'defensefinder-systems'], tsv]}.set{ ch_merge_systems }
    SYSTEMS_CONCAT(ch_merge_systems, 'tsv', 'tsv')
    ch_merged_systems = ch_merged_systems.mix(SYSTEMS_CONCAT.out.csv)
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
    versions = ch_versions // channel: [ versions.yml ]
}
