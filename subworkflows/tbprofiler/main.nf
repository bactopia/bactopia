//
// tbprofiler - Detect resistance and lineages of Mycobacterium tuberculosis genomes
//
nextflow.preview.types = true

include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'
include { flattenPaths       } from 'plugin/nf-bactopia'
include { gather             } from 'plugin/nf-bactopia'

workflow TBPROFILER {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    TBPROFILER_PROFILE(reads)
    TBPROFILER_COLLATE(gather(TBPROFILER_PROFILE.out.json, 'tbprofiler', 'json'))

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = TBPROFILER_PROFILE.out.csv
    json: Channel<Tuple<Map, Path>> = TBPROFILER_PROFILE.out.json
    txt: Channel<Tuple<Map, Path>> = TBPROFILER_PROFILE.out.txt
    bam: Channel<Tuple<Map, Path>> = TBPROFILER_PROFILE.out.bam
    vcf: Channel<Tuple<Map, Path>> = TBPROFILER_PROFILE.out.vcf
    merged_csv: Channel<Tuple<Map, Path>> = TBPROFILER_COLLATE.out.csv
    variants_csv: Channel<Tuple<Map, Path>> = TBPROFILER_COLLATE.out.variants_csv
    variants_txt: Channel<Tuple<Map, Path>> = TBPROFILER_COLLATE.out.variants_txt
    itol: Channel<Tuple<Map, Path>> = TBPROFILER_COLLATE.out.itol

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.csv,
        TBPROFILER_PROFILE.out.json,
        TBPROFILER_PROFILE.out.txt,
        TBPROFILER_PROFILE.out.bam,
        TBPROFILER_PROFILE.out.vcf,
        TBPROFILER_COLLATE.out.csv,
        TBPROFILER_COLLATE.out.variants_csv,
        TBPROFILER_COLLATE.out.variants_txt,
        TBPROFILER_COLLATE.out.itol
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.logs,
        TBPROFILER_COLLATE.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.nf_logs,
        TBPROFILER_COLLATE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.versions,
        TBPROFILER_COLLATE.out.versions
    ])
}
