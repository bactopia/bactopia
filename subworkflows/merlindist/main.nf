/**
 * Identify species from assembly and read data using Mash distances.
 *
 * This subworkflow performs rapid species identification using [Mash](https://github.com/marbl/Mash)
 * distance calculations against a reference database. It is a core component of the MERLIN
 * (MinER assisted species-specific bactopia tool seLectIoN) pipeline, responsible for determining
 * which species-specific typing tools should be run based on the detected organism. The workflow
 * outputs channels filtered by detected genera for downstream species-specific analysis.
 *
 * @status stable
 * @keywords species, identification, mash, distance, classification, taxonomy
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation mash
 *
 * @modules merlin_dist
 *
 * @input tuple(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format for species identification
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @input mash_db
 * Mash sketch database for rapid species identification
 *
 * @output dist                Mash distance results showing top species matches
 * @output escherichia         Samples identified as Escherichia/Shigella with assembly data
 * @output escherichia_fq      Samples identified as Escherichia/Shigella with read data (5-slot)
 * @output escherichia_fna_fq  Samples identified as Escherichia/Shigella with both assembly and read data (6-slot)
 * @output haemophilus         Samples identified as Haemophilus with assembly data
 * @output klebsiella          Samples identified as Klebsiella with assembly data
 * @output legionella          Samples identified as Legionella with assembly data
 * @output listeria            Samples identified as Listeria with assembly data
 * @output mycobacterium       Samples identified as Mycobacterium with assembly data
 * @output mycobacterium_fq    Samples identified as Mycobacterium with read data (5-slot)
 * @output neisseria           Samples identified as Neisseria with assembly data
 * @output pseudomonas         Samples identified as Pseudomonas with assembly data
 * @output salmonella          Samples identified as Salmonella with assembly data
 * @output salmonella_fq       Samples identified as Salmonella with read data (5-slot)
 * @output staphylococcus      Samples identified as Staphylococcus with assembly data
 * @output streptococcus       Samples identified as Streptococcus with assembly data
 * @output streptococcus_fq    Samples identified as Streptococcus with read data (5-slot)
 * @output results             Aggregated results channel containing all output files
 * @output logs                Aggregated logs channel containing all execution logs
 * @output nf_logs             Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions            Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MERLIN_DIST  } from '../../modules/merlin/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow MERLINDIST {
    take:
    ch_seqs: Channel<Tuple<Map, Path, Path?, Path?, Path?, Path?>>
    ch_mash_db: Path

    main:
    MERLIN_DIST(ch_seqs, ch_mash_db)

    emit:
    dist: Channel<Tuple<Map, Set<Path>>> = MERLIN_DIST.out.dist
    escherichia: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.escherichia
    escherichia_fq: Channel<Tuple<Map, Path?, Path?, Path?, Path?, Path?>> = MERLIN_DIST.out.escherichia_fq
    escherichia_fna_fq: Channel<Tuple<Map, Path, Path?, Path?, Path?, Path?, Path?>> = MERLIN_DIST.out.escherichia_fna_fq
    haemophilus: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.haemophilus
    klebsiella: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.klebsiella
    legionella: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.legionella
    listeria: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.listeria
    mycobacterium: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.mycobacterium
    mycobacterium_fq: Channel<Tuple<Map, Path?, Path?, Path?, Path?, Path?>> = MERLIN_DIST.out.mycobacterium_fq
    neisseria: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.neisseria
    pseudomonas: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.pseudomonas
    salmonella: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.salmonella
    salmonella_fq: Channel<Tuple<Map, Path?, Path?, Path?, Path?, Path?>> = MERLIN_DIST.out.salmonella_fq
    staphylococcus: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.staphylococcus
    streptococcus: Channel<Tuple<Map, Path, Path?>> = MERLIN_DIST.out.streptococcus
    streptococcus_fq: Channel<Tuple<Map, Path?, Path?, Path?, Path?, Path?>> = MERLIN_DIST.out.streptococcus_fq

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
