//
// mashdist - Calculate Mash distances between sequences
//
include { MASH_DIST } from '../../modules/mash/dist/main'
include { MERLIN_DIST } from '../../modules/mash/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MASHDIST {
    take:
    ch_seqs      // channel: [ val(meta), [ query ] ]
    ch_reference // channel: reference sketch file

    main:
    ch_versions = Channel.empty()
    ch_merged_mashdist = Channel.empty()

    // Run regular MASH_DIST (not Merlin mode)
    MASH_DIST(ch_seqs, ch_reference)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions.first())

    // Collect and concat results
    MASH_DIST.out.dist.collect{_meta, dist -> dist}.map{ dist -> [[id:'mashdist'], dist]}.set{ ch_merge_mashdist }
    CSVTK_CONCAT(ch_merge_mashdist, 'tsv', 'tsv')
    ch_merged_mashdist = ch_merged_mashdist.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    dist        = MASH_DIST.out.dist
    merged_dist = ch_merged_mashdist
    logs        = MASH_DIST.out.logs.mix(CSVTK_CONCAT.out.logs)
    nf_logs     = MASH_DIST.out.nf_begin.mix(MASH_DIST.out.nf_err, MASH_DIST.out.nf_log, MASH_DIST.out.nf_out, MASH_DIST.out.nf_run, MASH_DIST.out.nf_sh, MASH_DIST.out.nf_trace, CSVTK_CONCAT.out.nf_begin, CSVTK_CONCAT.out.nf_err, CSVTK_CONCAT.out.nf_log, CSVTK_CONCAT.out.nf_out, CSVTK_CONCAT.out.nf_run, CSVTK_CONCAT.out.nf_sh, CSVTK_CONCAT.out.nf_trace)
    versions    = ch_versions
}

workflow MERLINDIST {
    take:
    ch_seqs     // channel: [ val(meta), [ query, reads ] ]
    ch_mash_db  // channel: mash database

    main:
    ch_versions = Channel.empty()

    MERLIN_DIST(ch_seqs, ch_mash_db)
    ch_versions = ch_versions.mix(MERLIN_DIST.out.versions.first())

    emit:
    dist               = MERLIN_DIST.out.dist
    escherichia        = MERLIN_DIST.out.escherichia
    escherichia_fq     = MERLIN_DIST.out.escherichia_fq
    escherichia_fna_fq = MERLIN_DIST.out.escherichia_fna_fq
    haemophilus        = MERLIN_DIST.out.haemophilus
    klebsiella         = MERLIN_DIST.out.klebsiella
    legionella         = MERLIN_DIST.out.legionella
    listeria           = MERLIN_DIST.out.listeria
    mycobacterium      = MERLIN_DIST.out.mycobacterium
    mycobacterium_fq   = MERLIN_DIST.out.mycobacterium_fq
    neisseria          = MERLIN_DIST.out.neisseria
    pseudomonas        = MERLIN_DIST.out.pseudomonas
    salmonella         = MERLIN_DIST.out.salmonella
    salmonella_fq      = MERLIN_DIST.out.salmonella_fq
    staphylococcus     = MERLIN_DIST.out.staphylococcus
    streptococcus      = MERLIN_DIST.out.streptococcus
    streptococcus_fq   = MERLIN_DIST.out.streptococcus_fq
    logs               = MERLIN_DIST.out.logs
    nf_logs            = MERLIN_DIST.out.nf_begin.mix(MERLIN_DIST.out.nf_err, MERLIN_DIST.out.nf_log, MERLIN_DIST.out.nf_out, MERLIN_DIST.out.nf_run, MERLIN_DIST.out.nf_sh, MERLIN_DIST.out.nf_trace)
    versions           = ch_versions
}
