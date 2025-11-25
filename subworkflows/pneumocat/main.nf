//
// pneumocat - Assign capsular type to Streptococcus pneumoniae from sequence reads
//
nextflow.preview.types = true

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'

workflow PNEUMOCAT {
    take:
    fastq // channel: [ val(meta), [ fastq ] ]

    main:
    PNEUMOCAT_MODULE(fastq)

    emit:
    // Individual outputs
    xml = PNEUMOCAT_MODULE.out.xml
    txt = PNEUMOCAT_MODULE.out.txt

    // Generic aggregate outputs
    results = PNEUMOCAT_MODULE.out.xml.mix(
        PNEUMOCAT_MODULE.out.txt
    )
    logs = PNEUMOCAT_MODULE.out.logs
    nf_logs = PNEUMOCAT_MODULE.out.nf_begin.mix(
        PNEUMOCAT_MODULE.out.nf_err,
        PNEUMOCAT_MODULE.out.nf_log,
        PNEUMOCAT_MODULE.out.nf_out,
        PNEUMOCAT_MODULE.out.nf_run,
        PNEUMOCAT_MODULE.out.nf_sh,
        PNEUMOCAT_MODULE.out.nf_trace
    )
    versions = PNEUMOCAT_MODULE.out.versions
}
