//
// tbprofiler - Detect resistance and lineages of Mycobacterium tuberculosis genomes
//
include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'

workflow TBPROFILER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    TBPROFILER_PROFILE(reads)
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE.out.versions)
    ch_logs = ch_logs.mix(TBPROFILER_PROFILE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TBPROFILER_PROFILE.out.nf_logs)

    // Merge results
    TBPROFILER_PROFILE.out.json.collect{ _meta, json -> json }.map{ json -> [[id:'tbprofiler'], json] }.set{ ch_merge_tbprofiler }
    TBPROFILER_COLLATE(ch_merge_tbprofiler)
    ch_versions = ch_versions.mix(TBPROFILER_COLLATE.out.versions)
    ch_logs = ch_logs.mix(TBPROFILER_COLLATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TBPROFILER_COLLATE.out.nf_logs)

    emit:
    bam = TBPROFILER_PROFILE.out.bam
    csv = TBPROFILER_PROFILE.out.csv
    merged_csv = TBPROFILER_COLLATE.out.csv
    json = TBPROFILER_PROFILE.out.json
    txt = TBPROFILER_PROFILE.out.txt
    vcf = TBPROFILER_PROFILE.out.vcf
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
