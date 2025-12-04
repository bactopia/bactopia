//
// mashdist - Calculate Mash distances between sequences
//
nextflow.preview.types = true

include { MASH_DIST    } from '../../modules/mash/dist/main'
include { MERLIN_DIST  } from '../../modules/mash/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow MASHDIST {
    take:
    seqs: Channel<Tuple<Map, Set<Path>>>
    reference: Path

    main:
    MASH_DIST(seqs, reference)
    CSVTK_CONCAT(gather(MASH_DIST.out.dist, 'mashdist', 'dist'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    dist: Channel<Tuple<Map, Path>> = MASH_DIST.out.dist
    merged_dist: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.dist,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}

workflow MERLINDIST {
    take:
    ch_seqs: Channel<Tuple<Map, Set<Path>, Set<Path>>>
    ch_mash_db: Path

    main:
    MERLIN_DIST(ch_seqs, ch_mash_db)

    emit:
    dist: Channel<Tuple<Map, Set<Path>>> = MERLIN_DIST.out.dist
    escherichia: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.escherichia
    escherichia_fq: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.escherichia_fq
    escherichia_fna_fq: Channel<Tuple<Map, Set<Path>, Set<Path>, Path>> = MERLIN_DIST.out.escherichia_fna_fq
    haemophilus: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.haemophilus
    klebsiella: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.klebsiella
    legionella: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.legionella
    listeria: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.listeria
    mycobacterium: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.mycobacterium
    mycobacterium_fq: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.mycobacterium_fq
    neisseria: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.neisseria
    pseudomonas: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.pseudomonas
    salmonella: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.salmonella
    salmonella_fq: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.salmonella_fq
    staphylococcus: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.staphylococcus
    streptococcus: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.streptococcus
    streptococcus_fq: Channel<Tuple<Map, Set<Path>, Path>> = MERLIN_DIST.out.streptococcus_fq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MERLIN_DIST.out.dist,
        MERLIN_DIST.out.escherichia,
        MERLIN_DIST.out.escherichia_fq,
        MERLIN_DIST.out.escherichia_fna_fq,
        MERLIN_DIST.out.haemophilus,
        MERLIN_DIST.out.klebsiella,
        MERLIN_DIST.out.legionella,
        MERLIN_DIST.out.listeria,
        MERLIN_DIST.out.mycobacterium,
        MERLIN_DIST.out.mycobacterium_fq,
        MERLIN_DIST.out.neisseria,
        MERLIN_DIST.out.pseudomonas,
        MERLIN_DIST.out.salmonella,
        MERLIN_DIST.out.salmonella_fq,
        MERLIN_DIST.out.staphylococcus,
        MERLIN_DIST.out.streptococcus,
        MERLIN_DIST.out.streptococcus_fq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([MERLIN_DIST.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([MERLIN_DIST.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([MERLIN_DIST.out.versions])
}
