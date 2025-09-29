//
// qc - Quality control of Illumina and ONT reads
//
include { QC as QC_MODULE } from '../../modules/bactopia/qc/main'

workflow QC {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    adapters
    phix

    main:
    QC_MODULE(reads, adapters, phix)

    emit:
    // Individual outputs
    fastq = QC_MODULE.out.fastq
    fastq_only = QC_MODULE.out.fastq_only
    error = QC_MODULE.out.error

    // Generic aggregate outputs
    results = QC_MODULE.out.fastq_only.mix(
        QC_MODULE.out.error
    )
    logs = QC_MODULE.out.versions
    nf_logs = QC_MODULE.out.nf_begin.mix(
        QC_MODULE.out.nf_err,
        QC_MODULE.out.nf_log,
        QC_MODULE.out.nf_out,
        QC_MODULE.out.nf_run,
        QC_MODULE.out.nf_sh,
        QC_MODULE.out.nf_trace
    )
    versions = QC_MODULE.out.logs
}
