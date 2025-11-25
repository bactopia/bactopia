//
// scoary - GWAS analysis using pangenome outputs
//
nextflow.preview.types = true

include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'

workflow SCOARY {
    take:
    csv // channel: [ val(meta), [ csv ] ]
    traits

    main:    
    SCOARY_MODULE(csv, traits)

    emit:
    // Individual outputs
    csv = SCOARY_MODULE.out.csv

    // Generic aggregate outputs
    results = SCOARY_MODULE.out.csv
    logs = SCOARY_MODULE.out.logs
    nf_logs = SCOARY_MODULE.out.nf_begin.mix(
        SCOARY_MODULE.out.nf_err,
        SCOARY_MODULE.out.nf_log,
        SCOARY_MODULE.out.nf_out,
        SCOARY_MODULE.out.nf_run,
        SCOARY_MODULE.out.nf_sh,
        SCOARY_MODULE.out.nf_trace
    )
    versions = SCOARY_MODULE.out.versions
}
