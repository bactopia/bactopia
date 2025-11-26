//
// ncbigenomedownload - Quickly download assemblies from NCBI's Assembly database
//
nextflow.preview.types = true

include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../modules/ncbigenomedownload/main'
include { flattenPaths                                    } from 'plugin/nf-bactopia'
include { gather                                          } from 'plugin/nf-bactopia'

workflow NCBIGENOMEDOWNLOAD {

    take:
    accessions: Channel<String>

    main:
    NCBIGENOMEDOWNLOAD_MODULE(accessions)
    ch_to_bactopia_tools = NCBIGENOMEDOWNLOAD_MODULE.out.all.map { path -> [[id: file(path).getSimpleName()], file(path)] }

    emit:
    // Individual outputs
    bactopia_tools: Channel<Tuple<Map, Path>> = ch_to_bactopia_tools
    gbk: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gbk
    fna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.fna
    rm: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rm
    features: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.features
    gff: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gff
    faa: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.faa
    gpff: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gpff
    wgs_gbk: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk
    cds: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.cds
    rna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rna
    rna_fna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna
    report: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.report
    stats: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.stats

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        NCBIGENOMEDOWNLOAD_MODULE.out.gbk,
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
    ])
    logs: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.versions
}
