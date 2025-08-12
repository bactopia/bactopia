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
    
    MOBSUITE_RECON(fasta)
    ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions.first())
    ch_logs = ch_logs.mix(MOBSUITE_RECON.out.logs)
    
    // Merge results
    MOBSUITE_RECON.out.mobtyper_results.collect{_meta, summary -> summary}.map{ summary -> [[id:'mobsuite'], summary]}.set{ ch_merge_mobsuite }
    CSVTK_CONCAT(ch_merge_mobsuite, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    merged_reports = CSVTK_CONCAT.out.csv
    chromosome = MOBSUITE_RECON.out.chromosome
    contig_report = MOBSUITE_RECON.out.contig_report
    mobtyper_results = MOBSUITE_RECON.out.mobtyper_results
    plasmids = MOBSUITE_RECON.out.plasmids
    logs = ch_logs
    nf_logs = MOBSUITE_RECON.out.nf_begin.mix(
        MOBSUITE_RECON.out.nf_err,
        MOBSUITE_RECON.out.nf_log,
        MOBSUITE_RECON.out.nf_out,
        MOBSUITE_RECON.out.nf_run,
        MOBSUITE_RECON.out.nf_sh,
        MOBSUITE_RECON.out.nf_trace,
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
