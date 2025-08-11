//
// pirate - Pangenome toolbox for bacterial genomes
//
include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'

workflow PIRATE {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    PIRATE_MODULE(gff)
    ch_versions = ch_versions.mix(PIRATE_MODULE.out.versions)

    emit:
    aln = PIRATE_MODULE.out.aln
    csv = PIRATE_MODULE.out.csv
    logs = PIRATE_MODULE.out.logs
    nf_logs = PIRATE_MODULE.out.nf_begin.mix(
        PIRATE_MODULE.out.nf_err,
        PIRATE_MODULE.out.nf_log,
        PIRATE_MODULE.out.nf_out,
        PIRATE_MODULE.out.nf_run,
        PIRATE_MODULE.out.nf_sh,
        PIRATE_MODULE.out.nf_trace
    )
    versions = ch_versions
}
