//
// mashdist - Calculate Mash distances between sequences
//
nextflow.preview.types = true

include { MASH_DIST } from '../../modules/mash/dist/main'
include { MERLIN_DIST } from '../../modules/mash/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MASHDIST {
    take:
    seqs      // channel: [ val(meta), [ query ] ]
    reference // channel: reference sketch file

    main:
    MASH_DIST(seqs, reference)

    // Merge results
    ch_merge_mashdist = MASH_DIST.out.dist.collect{_meta, dist -> dist}.map{ dist -> [[id:'mashdist'], dist]}
    CSVTK_CONCAT(ch_merge_mashdist, 'tsv', 'tsv')

    emit:
    // Individual outputs
    dist = MASH_DIST.out.dist
    merged_dist = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results  = MASH_DIST.out.dist.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = MASH_DIST.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs  = MASH_DIST.out.nf_begin.mix(
        MASH_DIST.out.nf_err,
        MASH_DIST.out.nf_log,
        MASH_DIST.out.nf_out,
        MASH_DIST.out.nf_run,
        MASH_DIST.out.nf_sh,
        MASH_DIST.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = MASH_DIST.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}

workflow MERLINDIST {
    take:
    ch_seqs     // channel: [ val(meta), [ query, reads ] ]
    ch_mash_db  // channel: mash database

    main:
    MERLIN_DIST(ch_seqs, ch_mash_db)

    emit:
    // Individual outputs
    dist = MERLIN_DIST.out.dist
    escherichia = MERLIN_DIST.out.escherichia
    escherichia_fq = MERLIN_DIST.out.escherichia_fq
    escherichia_fna_fq = MERLIN_DIST.out.escherichia_fna_fq
    haemophilus = MERLIN_DIST.out.haemophilus
    klebsiella = MERLIN_DIST.out.klebsiella
    legionella = MERLIN_DIST.out.legionella
    listeria = MERLIN_DIST.out.listeria
    mycobacterium = MERLIN_DIST.out.mycobacterium
    mycobacterium_fq = MERLIN_DIST.out.mycobacterium_fq
    neisseria = MERLIN_DIST.out.neisseria
    pseudomonas = MERLIN_DIST.out.pseudomonas
    salmonella = MERLIN_DIST.out.salmonella
    salmonella_fq = MERLIN_DIST.out.salmonella_fq
    staphylococcus = MERLIN_DIST.out.staphylococcus
    streptococcus = MERLIN_DIST.out.streptococcus
    streptococcus_fq = MERLIN_DIST.out.streptococcus_fq

    // Generic aggregate outputs
    results = MERLIN_DIST.out.dist
    logs = MERLIN_DIST.out.logs
    nf_logs = MERLIN_DIST.out.nf_begin.mix(
        MERLIN_DIST.out.nf_err,
        MERLIN_DIST.out.nf_log,
        MERLIN_DIST.out.nf_out,
        MERLIN_DIST.out.nf_run,
        MERLIN_DIST.out.nf_sh,
        MERLIN_DIST.out.nf_trace
    )
    versions = MERLIN_DIST.out.versions
}
