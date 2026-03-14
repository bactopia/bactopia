/**
 * Identify transposase insertion sites in bacterial genomes.
 *
 * This subworkflow maps insertion sequence (IS) positions in bacterial genomes
 * using [ISMapper](https://github.com/jhawkey/IS_mapper). The tool identifies
 * transposase insertion sites from short read sequence data by mapping reads
 * to reference sequences and detecting insertion sites with high precision.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords insertion, sequence, transposase, mobile genetic elements
 * @tags complexity:simple input-type:multiple output-type:single
 * @citation ismapper
 *
 * @modules ismapper
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by ISMapper)
 * - `lr`: Long reads (not supported by ISMapper)
 *
 * @input reference
 * Reference genome in FASTA format for mapping
 *
 * @input insertions
 * Insertion sequence reference file containing IS elements to map
 *
 * @output sample_outputs
 * - `supplemental`: Directory containing the final tables of insertion sites and visual summaries
 */
nextflow.preview.types = true

include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'

workflow ISMAPPER {
    take:
    reads: Channel<Record>
    reference: Path
    insertions: Path

    main:
    ISMAPPER_MODULE(reads, reference, insertions)

    emit:
    sample_outputs = ISMAPPER_MODULE.out
}
