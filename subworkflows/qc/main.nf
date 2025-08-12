//
// qc - Quality control of Illumina and ONT reads
//
include { QC as QC_MODULE } from '../../modules/bactopia/qc/main'

workflow QC {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    // Set up input files


    // QC reads
    QC_MODULE(reads, ADAPTERS, PHIX)
    ch_versions = ch_versions.mix(QC_MODULE.out.versions)
    ch_logs = ch_logs.mix(QC_MODULE.out.logs)

    emit:

    emit:
    ADAPTERS = params.adapters ? file(params.adapters) : []
    PHIX = params.phix ? file(params.phix) : []
    fastq = QC_MODULE.out.fastq
    fastq_only = QC_MODULE.out.fastq_only
    logs = ch_logs
    nf_logs = QC_MODULE.out.nf_begin.mix(
        QC_MODULE.out.nf_err,
        QC_MODULE.out.nf_log,
        QC_MODULE.out.nf_out,
        QC_MODULE.out.nf_run,
        QC_MODULE.out.nf_sh,
        QC_MODULE.out.nf_trace
    )
    versions = ch_versions
}
