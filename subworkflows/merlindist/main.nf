/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat, merlin_dist, mash_dist
 *
 * @input seqs
 * Channel containing seqs data
 *
 * @input reference
 * Channel containing reference data
 *
 * @output dist        Dist
 * @output merged_dist Merged Dist
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MERLIN_DIST  } from '../../modules/merlin/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

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
