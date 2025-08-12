//
// pneumocat - Assign capsular type to Streptococcus pneumoniae from sequence reads
//
include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'

workflow PNEUMOCAT {
    take:
    fastq // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    PNEUMOCAT_MODULE(fastq)
    ch_versions = ch_versions.mix(PNEUMOCAT_MODULE.out.versions.first()

    emit:
    xml = PNEUMOCAT_MODULE.out.xml
    logs = PNEUMOCAT_MODULE.out.logs
    nf_logs = PNEUMOCAT_MODULE.out.nf_begin.mix(
        PNEUMOCAT_MODULE.out.nf_err,
        PNEUMOCAT_MODULE.out.nf_log,
        PNEUMOCAT_MODULE.out.nf_out,
        PNEUMOCAT_MODULE.out.nf_run,
        PNEUMOCAT_MODULE.out.nf_sh,
        PNEUMOCAT_MODULE.out.nf_trace
    )
    versions = ch_versions
}
