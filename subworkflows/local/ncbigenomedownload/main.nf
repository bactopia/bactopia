//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
download_opts = [
    params.kingdom,
    "--section ${params.section}",
    "--formats ${params.format}",
    "--assembly-levels ${params.assembly_level}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../../modules/nf-core/modules/ncbigenomedownload/main' addParams( options: [ args: "${download_opts}", is_module: true] )

workflow NCBIGENOMEDOWNLOAD {
    main:
    ch_versions = Channel.empty()

    inputs = [[
        limit: params.limit,
        has_accessions: (params.accessions ? true, false),
        accession: params.accession,
        species: params.species
    ]]
    accessions = params.accessions ? file(params.accessions), []
    NCBIGENOMEDOWNLOAD_MODULE(inputs, accessions)
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD_MODULE.out.versions.first())

    emit:
    NCBIGENOMEDOWNLOAD_MODULE.out.gbk
    NCBIGENOMEDOWNLOAD_MODULE.out.fna
    NCBIGENOMEDOWNLOAD_MODULE.out.rm
    NCBIGENOMEDOWNLOAD_MODULE.out.features
    NCBIGENOMEDOWNLOAD_MODULE.out.gff
    NCBIGENOMEDOWNLOAD_MODULE.out.faa
    NCBIGENOMEDOWNLOAD_MODULE.out.gpff
    NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk
    NCBIGENOMEDOWNLOAD_MODULE.out.cds
    NCBIGENOMEDOWNLOAD_MODULE.out.rna
    NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna
    NCBIGENOMEDOWNLOAD_MODULE.out.report
    NCBIGENOMEDOWNLOAD_MODULE.out.stats
    versions = ch_versions // channel: [ versions.yml ]
}
