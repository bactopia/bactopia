//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../modules/ncbigenomedownload/main'

workflow NCBIGENOMEDOWNLOAD {
    main:
    ch_versions = Channel.empty()
    
    META = [
        id: "",
        limit: params.limit,
        has_accessions: params.accessions ? true : false,
        accession: params.accession,
        species: params.species
    ]
    ACCESSIONS = params.accessions ? file(params.accessions) : []

    NCBIGENOMEDOWNLOAD_MODULE(META, ACCESSIONS)
    NCBIGENOMEDOWNLOAD_MODULE.out.all.flatten().map{ [[id: file(it).getSimpleName()], file(it)]}.set{ ch_to_bactopia_tools }
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD_MODULE.out.versions)

    emit:
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
    logs = NCBIGENOMEDOWNLOAD_MODULE.out.logs
    nf_logs = NCBIGENOMEDOWNLOAD_MODULE.out.nf_begin.mix(
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_err,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_log,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_out,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_run,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_sh,
        NCBIGENOMEDOWNLOAD_MODULE.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
