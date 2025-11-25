//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
nextflow.preview.types = true

include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../modules/ncbigenomedownload/main'

workflow NCBIGENOMEDOWNLOAD {

    take:
    accessions

    main:
    NCBIGENOMEDOWNLOAD_MODULE(accessions)
    ch_to_bactopia_tools = NCBIGENOMEDOWNLOAD_MODULE.out.all.map { path -> [[id: file(path).getSimpleName()], file(path)] }

    emit:
    // Individual outputs
    bactopia_tools = ch_to_bactopia_tools
    gbk = NCBIGENOMEDOWNLOAD_MODULE.out.gbk
    fna = NCBIGENOMEDOWNLOAD_MODULE.out.fna
    rm = NCBIGENOMEDOWNLOAD_MODULE.out.rm
    features = NCBIGENOMEDOWNLOAD_MODULE.out.features
    gff = NCBIGENOMEDOWNLOAD_MODULE.out.gff
    faa = NCBIGENOMEDOWNLOAD_MODULE.out.faa
    gpff = NCBIGENOMEDOWNLOAD_MODULE.out.gpff
    wgs_gbk = NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk
    cds = NCBIGENOMEDOWNLOAD_MODULE.out.cds
    rna = NCBIGENOMEDOWNLOAD_MODULE.out.rna
    rna_fna = NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna
    report = NCBIGENOMEDOWNLOAD_MODULE.out.report
    stats = NCBIGENOMEDOWNLOAD_MODULE.out.stats

    // Generic aggregate outputs
    results = NCBIGENOMEDOWNLOAD_MODULE.out.gbk.mix(
        NCBIGENOMEDOWNLOAD_MODULE.out.fna,
        NCBIGENOMEDOWNLOAD_MODULE.out.rm,
        NCBIGENOMEDOWNLOAD_MODULE.out.features,
        NCBIGENOMEDOWNLOAD_MODULE.out.gff,
        NCBIGENOMEDOWNLOAD_MODULE.out.faa,
        NCBIGENOMEDOWNLOAD_MODULE.out.gpff,
        NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk,
        NCBIGENOMEDOWNLOAD_MODULE.out.cds,
        NCBIGENOMEDOWNLOAD_MODULE.out.rna,
        NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna,
        NCBIGENOMEDOWNLOAD_MODULE.out.report,
        NCBIGENOMEDOWNLOAD_MODULE.out.stats
    )
    logs = NCBIGENOMEDOWNLOAD_MODULE.out.logs
    nf_logs = NCBIGENOMEDOWNLOAD_MODULE.out.nf_begin.mix(
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_err,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_log,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_out,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_run,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_sh,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_trace
    )
    versions = NCBIGENOMEDOWNLOAD_MODULE.out.versions
}
