/**
 * Create genomic sketches and perform rapid taxonomic classification.
 *
 * This subworkflow generates MinHash sketches from assembled genomes using [Mash](https://github.com/marbl/Mash)
 * and [Sourmash](https://github.com/dib-lab/sourmash). The sketches are compared against reference databases
 * to identify taxonomic classification and find closely related genomes. Mash queries against RefSeq while
 * Sourmash uses the GTDB database for comprehensive taxonomic placement.
 *
 * @status stable
 * @keywords taxonomy, classification, minhash, sketch, mash, sourmash, refseq, gtdb
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, compression
 * @citation mash, sourmash
 *
 * @modules sketcher
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input mash_db
 * Path to the Mash RefSeq database for taxonomic classification
 *
 * @input sourmash_db
 * Path to the Sourmash GTDB LCA database for taxonomic classification
 *
 * @output sig       Sourmash signature files for downstream comparative analyses
 * @output msh       Mash sketch files created at k=21 and k=31
 * @output mash      Classification report from Mash Screen against RefSeq database
 * @output sourmash  Classification report from Sourmash LCA against GTDB database
 * @output results   Aggregated results channel containing all output files
 * @output logs      Aggregated logs channel containing all execution logs
 * @output nf_logs   Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions  Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SKETCHER as SKETCHER_MODULE } from '../../../modules/bactopia/sketcher/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SKETCHER {
    take:
    assembly: Channel<Tuple<Map, Path>>
    mash_db: Path
    sourmash_db: Path

    main:
    SKETCHER_MODULE(assembly, mash_db, sourmash_db)

    emit:
    // Individual outputs
    sig: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sig
    msh: Channel<Tuple<Map, Set<Path>>> = SKETCHER_MODULE.out.msh
    mash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.mash
    sourmash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sourmash

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SKETCHER_MODULE.out.sig,
        SKETCHER_MODULE.out.msh,
        SKETCHER_MODULE.out.mash,
        SKETCHER_MODULE.out.sourmash
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.versions])
}
