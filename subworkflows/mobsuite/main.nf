//
// mobsuite - Reconstruct and annotate plasmids in bacterial assemblies
//
include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MOBSUITE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    MOBSUITE_RECON(fasta)
    ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions.first())
    ch_logs = ch_logs.mix(MOBSUITE_RECON.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MOBSUITE_RECON.out.nf_logs)

    // Merge results
    MOBSUITE_RECON.out.mobtyper_results.collect{ _meta, summary -> summary }.map{ summary -> [[id:'mobsuite'], summary] }.set{ ch_merge_mobsuite }
    CSVTK_CONCAT(ch_merge_mobsuite, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    chromosome = MOBSUITE_RECON.out.chromosome
    contig_report = MOBSUITE_RECON.out.contig_report
    plasmids = MOBSUITE_RECON.out.plasmids
    mobtyper_results = MOBSUITE_RECON.out.mobtyper_results
    merged_reports = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
