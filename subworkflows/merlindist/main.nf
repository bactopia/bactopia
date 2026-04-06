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
 * @input record(meta, fna, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format for species identification
 * - `r1?`: Illumina R1 reads (paired-end) or null
 * - `r2?`: Illumina R2 reads (paired-end) or null
 * - `se?`: Single-end Illumina reads or null
 * - `lr?`: Long reads (ONT/PacBio) or null
 *
 * @input mash_db
 * Mash sketch database for rapid species identification
 *
 * @output sample_outputs
 * - `dist`: The raw Mash distance results
 * - `fna`: Passthrough of assembled contigs
 * - `r1`: Passthrough of Illumina R1 reads
 * - `r2`: Passthrough of Illumina R2 reads
 * - `se`: Passthrough of single-end reads
 * - `lr`: Passthrough of long reads
 * - `escherichia`: Conditional marker file triggering Escherichia analysis tools
 * - `haemophilus`: Conditional marker file triggering Haemophilus analysis tools
 * - `klebsiella`: Conditional marker file triggering Klebsiella analysis tools
 * - `legionella`: Conditional marker file triggering Legionella analysis tools
 * - `listeria`: Conditional marker file triggering Listeria analysis tools
 * - `mycobacterium`: Conditional marker file triggering Mycobacterium analysis tools
 * - `neisseria`: Conditional marker file triggering Neisseria analysis tools
 * - `pseudomonas`: Conditional marker file triggering Pseudomonas analysis tools
 * - `salmonella`: Conditional marker file triggering Salmonella analysis tools
 * - `staphylococcus`: Conditional marker file triggering Staphylococcus analysis tools
 * - `streptococcus`: Conditional marker file triggering Streptococcus analysis tools
 * - `genus`: A marker file indicating the detected genus (for debugging)
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { MERLIN_DIST } from '../../modules/merlin/dist/main'
include { gatherCsvtk } from 'plugin/nf-bactopia'

workflow MERLINDIST {
    take:
    ch_seqs: Channel<Record>
    ch_mash_db: Value<Path>

    main:
    MERLIN_DIST(ch_seqs, ch_mash_db)

    emit:
    // Published outputs
    sample_outputs = MERLIN_DIST.out
    run_outputs = channel.empty()
}
