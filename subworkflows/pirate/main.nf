//
// pirate - Pangenome toolbox for bacterial genomes
//
nextflow.preview.types = true

include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'

workflow PIRATE {
    take:
    gff: Channel<Tuple<Map, Path>>

    main:
    PIRATE_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.aln
    csv: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.supplemental.mix(
        PIRATE_MODULE.out.aln,
        PIRATE_MODULE.out.csv
    )
    logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.nf_begin.mix(
        PIRATE_MODULE.out.nf_err,
        PIRATE_MODULE.out.nf_log,
        PIRATE_MODULE.out.nf_out,
        PIRATE_MODULE.out.nf_run,
        PIRATE_MODULE.out.nf_sh,
        PIRATE_MODULE.out.nf_trace
    )
    versions: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.versions
}
