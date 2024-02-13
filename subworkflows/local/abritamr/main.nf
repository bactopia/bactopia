//
// abritamr - A NATA accredited tool for reporting the presence of antimicrobial resistance genes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'abritamr')
options.args = [
    params.abritamr_identity ? "--identity ${params.abritamr_identity} " : "",
    params.abritamr_species ? "--species ${params.abritamr_species} " : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ABRITAMR_RUN } from '../../../modules/nf-core/abritamr/run/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'abritamr-concat', process_name: params.merge_folder] )

workflow ABRITAMR {
    take:
    fasta // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_summaries = Channel.empty()

    // Run AMRFinder
    ABRITAMR_RUN ( fasta )
    ch_versions = ch_versions.mix(ABRITAMR_RUN.out.versions.first())

    // Merge results
    ABRITAMR_RUN.out.summary.collect{meta, summary -> summary}.map{ summary -> [[id:'abritamr'], summary]}.set{ ch_merge_summary }
    CSVTK_CONCAT(ch_merge_summary, 'tsv', 'tsv')
    ch_merged_summaries = ch_merged_summaries.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    summary_tsv = ABRITAMR_RUN.out.summary
    merged_summary_tsv = ch_merged_summaries
    matches_tsv = ABRITAMR_RUN.out.matches
    partials_tsv = ABRITAMR_RUN.out.partials
    virulence_tsv = ABRITAMR_RUN.out.virulence
    amrfinder_tsv = ABRITAMR_RUN.out.amrfinder
    versions = ch_versions // channel: [ versions.yml ]
}
