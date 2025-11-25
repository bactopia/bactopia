//
// snpdists - Pairwise SNP distance matrix from a FASTA sequence alignment
//
nextflow.preview.types = true

include { SNPDISTS as SNPDISTS_MODULE } from '../../modules/snpdists/main'

workflow SNPDISTS {
    take:
    alignment: Channel<Tuple<Map, Path>>

    main:
    SNPDISTS_MODULE(alignment)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.tsv
    logs: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.nf_begin.mix(
        SNPDISTS_MODULE.out.nf_err,
        SNPDISTS_MODULE.out.nf_log,
        SNPDISTS_MODULE.out.nf_out,
        SNPDISTS_MODULE.out.nf_run,
        SNPDISTS_MODULE.out.nf_sh,
        SNPDISTS_MODULE.out.nf_trace
    )
    versions: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.versions
}
