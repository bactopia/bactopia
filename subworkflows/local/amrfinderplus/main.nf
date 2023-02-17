//
// amrfinderplus - Identify antimicrobial resistance in genes or proteins
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus')
options.args = [
    params.report_common ? "--report_common" : "",
    params.report_all_equal ? "--report_all_equal" : "",
    params.organism ? "--organism ${params.organism}" : "",
    "--ident_min ${params.ident_min}",
    "--coverage_min ${params.coverage_min}",
    "--translation_table ${params.translation_table}",
    "${params.amrfinder_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { AMRFINDERPLUS_RUN } from '../../../modules/nf-core/amrfinderplus/run/main' addParams( options: options )
include { CSVTK_CONCAT as GENES_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'amrfinderplus-genes-concat', process_name: params.merge_folder] )
include { CSVTK_CONCAT as PROTEINS_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'amrfinderplus-proteins-concat', process_name: params.merge_folder] )

workflow AMRFINDERPLUS {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ amrfinderplus_db ]

    main:
    ch_versions = Channel.empty()
    ch_amrfinder_db = Channel.empty()
    ch_merged_gene_reports = Channel.empty()
    ch_merged_protein_reports = Channel.empty()

    // Run AMRFinder
    AMRFINDERPLUS_RUN ( fasta, db )
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())

    // Merge results
    AMRFINDERPLUS_RUN.out.gene_report.collect{meta, report -> report}.map{ report -> [[id:'amrfinderplus-genes'], report]}.set{ ch_merge_gene_report }
    GENES_CONCAT(ch_merge_gene_report, 'tsv', 'tsv')
    ch_merged_gene_reports = ch_merged_gene_reports.mix(GENES_CONCAT.out.csv)
    ch_versions = ch_versions.mix(GENES_CONCAT.out.versions)

    AMRFINDERPLUS_RUN.out.protein_report.collect{meta, report -> report}.map{ report -> [[id:'amrfinderplus-proteins'], report]}.set{ ch_merge_protein_report }
    PROTEINS_CONCAT(ch_merge_protein_report, 'tsv', 'tsv')
    ch_merged_protein_reports = ch_merged_protein_reports.mix(PROTEINS_CONCAT.out.csv)
    ch_versions = ch_versions.mix(PROTEINS_CONCAT.out.versions)

    emit:
    gene_tsv = AMRFINDERPLUS_RUN.out.gene_report
    merged_gene_tsv = ch_merged_gene_reports
    protein_tsv = AMRFINDERPLUS_RUN.out.protein_report
    merged_protein_tsv = ch_merged_protein_reports
    mutation_reports = AMRFINDERPLUS_RUN.out.mutation_reports
    db = ch_amrfinder_db
    versions = ch_versions // channel: [ versions.yml ]
}
