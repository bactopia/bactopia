//
// tbprofiler - Detect resistance and lineages of Mycobacterium tuberculosis genomes
//
nextflow.preview.types = true

include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'

workflow TBPROFILER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    TBPROFILER_PROFILE(reads)

    // Merge results
    ch_merge_tbprofiler = TBPROFILER_PROFILE.out.json.collect{_meta, json -> json}.map{ json -> [[id:'tbprofiler'], json]}
    TBPROFILER_COLLATE(ch_merge_tbprofiler)

    emit:
    // Individual outputs
    csv = TBPROFILER_PROFILE.out.csv
    json = TBPROFILER_PROFILE.out.json
    txt = TBPROFILER_PROFILE.out.txt
    bam = TBPROFILER_PROFILE.out.bam
    vcf = TBPROFILER_PROFILE.out.vcf
    merged_csv = TBPROFILER_COLLATE.out.csv
    variants_csv = TBPROFILER_COLLATE.out.variants_csv
    variants_txt = TBPROFILER_COLLATE.out.variants_txt
    itol = TBPROFILER_COLLATE.out.itol

    // Generic aggregate outputs
    results = TBPROFILER_PROFILE.out.csv.mix(
        TBPROFILER_PROFILE.out.json,
        TBPROFILER_PROFILE.out.txt,
        TBPROFILER_PROFILE.out.bam,
        TBPROFILER_PROFILE.out.vcf,
        TBPROFILER_COLLATE.out.csv,
        TBPROFILER_COLLATE.out.variants_csv,
        TBPROFILER_COLLATE.out.variants_txt,
        TBPROFILER_COLLATE.out.itol
    )
    logs = TBPROFILER_PROFILE.out.logs.mix(
        TBPROFILER_COLLATE.out.logs
    )
    nf_logs = TBPROFILER_PROFILE.out.nf_begin.mix(
        TBPROFILER_PROFILE.out.nf_err,
        TBPROFILER_PROFILE.out.nf_log,
        TBPROFILER_PROFILE.out.nf_out,
        TBPROFILER_PROFILE.out.nf_run,
        TBPROFILER_PROFILE.out.nf_sh,
        TBPROFILER_PROFILE.out.nf_trace
    )
    versions = TBPROFILER_PROFILE.out.versions.mix(
        TBPROFILER_COLLATE.out.versions
    )
}
