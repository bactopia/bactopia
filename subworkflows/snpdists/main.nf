//
// snpdists - Pairwise SNP distance matrix from a FASTA sequence alignment
//
include { SNPDISTS as SNPDISTS_MODULE } from '../../modules/snpdists/main'

workflow SNPDISTS {
    take:
    alignment // channel: [ val(meta), [ alignment ] ]

    main:
    SNPDISTS_MODULE(alignment)

    emit:
    // Individual outputs
    tsv = SNPDISTS_MODULE.out.tsv

    // Generic aggregate outputs
    results = SNPDISTS_MODULE.out.tsv
    logs = SNPDISTS_MODULE.out.logs
    nf_logs = SNPDISTS_MODULE.out.nf_begin.mix(
        SNPDISTS_MODULE.out.nf_err,
        SNPDISTS_MODULE.out.nf_log,
        SNPDISTS_MODULE.out.nf_out,
        SNPDISTS_MODULE.out.nf_run,
        SNPDISTS_MODULE.out.nf_sh,
        SNPDISTS_MODULE.out.nf_trace
    )
    versions = SNPDISTS_MODULE.out.versions
}
