//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ncbigenomedownload')
options.is_module = params.wf == 'ncbigenomedownload' ? true : false
options.args = [
    params.kingdom,
    "--section ${params.section}",
    "--formats ${params.format}",
    "--assembly-levels ${params.assembly_level}",
    "--verbose",
    params.enable_conda ? "" : "--no-cache" 
].join(' ').replaceAll("\\s{2,}", " ").trim()
META = [
    id: "",
    limit: params.limit,
    has_accessions: params.accessions ? true : false,
    accession: params.accession,
    species: params.species
]
ACCESSIONS = params.accessions ? file(params.accessions) : []
include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../../modules/nf-core/ncbigenomedownload/main' addParams( options: options )

workflow NCBIGENOMEDOWNLOAD {
    main:
    ch_versions = Channel.empty()

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
    versions = ch_versions // channel: [ versions.yml ]
}
