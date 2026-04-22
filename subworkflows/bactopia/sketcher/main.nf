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
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation mash, sourmash
 *
 * @modules bactopia_sketcher
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input mash_db
 * Path to the Mash RefSeq database for taxonomic classification
 *
 * @input sourmash_db
 * Path to the Sourmash GTDB LCA database for taxonomic classification
 *
 * @output sample_outputs
 *   - `sig`: Sourmash signature file
 *   - `msh`: Mash sketch files for k=21 and k=31
 *   - `mash`: Mash Screen classification report against RefSeq
 *   - `sourmash`: Sourmash LCA classification report against GTDB
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { SKETCHER as SKETCHER_MODULE } from '../../../modules/bactopia/sketcher/main'

workflow SKETCHER {
    take:
    assembly: Channel<Record>
    mash_db: Path
    sourmash_db: Path

    main:
    ch_sketcher = SKETCHER_MODULE(assembly, mash_db, sourmash_db)

    emit:
    // Published outputs
    sample_outputs = ch_sketcher
    run_outputs = channel.empty()
}
