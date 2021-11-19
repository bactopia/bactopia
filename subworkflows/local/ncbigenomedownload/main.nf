//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
download_opts = [
    params.kingdom,
    "--section ${params.section}",
    "--formats ${params.format}",
    "--assembly-levels ${params.assembly_level}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
META = [
    id: "",
    limit: params.limit,
    has_accessions: params.accessions ? true : false,
    accession: params.accession,
    species: params.species
]
ACCESSIONS = params.accessions ? file(params.accessions) : []
include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../../modules/nf-core/modules/ncbigenomedownload/main' addParams( options: [ args: "${download_opts}", is_module: true] )

workflow NCBIGENOMEDOWNLOAD {
    main:
    ch_versions = Channel.empty()

    NCBIGENOMEDOWNLOAD_MODULE(META, ACCESSIONS)
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD_MODULE.out.versions.first())

    emit:
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
    versions = ch_versions // channel: [ versions.yml ]
}
