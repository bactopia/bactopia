//
// mobsuite - Reconstruct and annotate plasmids in bacterial assemblies
//
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MOBSUITE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    MOBSUITE_RECON(fasta)
    
    // Merge results
    ch_merge_mobsuite = MOBSUITE_RECON.out.txt.collect{_meta, summary -> summary}.map{ summary -> [[id:'mobsuite'], summary]}
    CSVTK_CONCAT(ch_merge_mobsuite, 'tsv', 'tsv')

    emit:
    txt = MOBSUITE_RECON.out.txt
    merged_reports = CSVTK_CONCAT.out.csv
    chromosome = MOBSUITE_RECON.out.chromosome
    contig_report = MOBSUITE_RECON.out.contig_report
    plasmids = MOBSUITE_RECON.out.plasmids

    // Generic aggregate outputs
    results = MOBSUITE_RECON.out.txt.mix(
        CSVTK_CONCAT.out.csv,
        MOBSUITE_RECON.out.chromosome,
        MOBSUITE_RECON.out.contig_report,
        MOBSUITE_RECON.out.plasmids
    )
    logs = MOBSUITE_RECON.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = MOBSUITE_RECON.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
