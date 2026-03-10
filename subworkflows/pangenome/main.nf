/**
 * Perform pangenome analysis with optional core-genome phylogeny.
 *
 * This subworkflow creates a pangenome from GFF3 annotation files using one of three
 * tools: [Panaroo](https://github.com/gtonkinhill/panaroo) (default),
 * [PIRATE](https://github.com/SionBayliss/PIRATE), or
 * [Roary](https://github.com/sanger-pathogens/roary). It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations using
 * [snp-dists](https://github.com/tseemann/snp-dists). The workflow conditionally executes
 * the selected pangenome tool based on Boolean parameters.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, conditional-logic
 * @citation pirate, panaroo, roary, snpdists
 *
 * @subworkflows pirate, roary, panaroo, snpdists
 *
 * @input tuple(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: Set of GFF3 annotation files from assembled genomes
 *
 * @input use_pirate
 * Boolean flag to use PIRATE for pangenome analysis
 *
 * @input use_roary
 * Boolean flag to use Roary for pangenome analysis
 *
 * @output sample_outputs
 * - `aln`: Core-genome alignment in FASTA format
 * - `csv`: Gene presence/absence matrix
 * - `supplemental`: Intermediate files and detailed outputs
 *
 * @output snpdists_outputs
 * - `tsv`: Pairwise SNP distance matrix from core-genome alignment
 */
nextflow.preview.types = true

include { PIRATE   } from '../pirate/main'
include { ROARY    } from '../roary/main'
include { PANAROO  } from '../panaroo/main'
include { SNPDISTS } from '../snpdists/main'

workflow PANGENOME {
    take:
    gff        : Channel<Record>
    use_pirate : Boolean
    use_roary  : Boolean

    main:
    ch_sample_outputs = channel.empty()
    ch_aln = channel.empty()

    // Choose pangenome tool based on params
    if (use_pirate) {
        PIRATE(gff)
        ch_sample_outputs = PIRATE.out.sample_outputs
        ch_aln = PIRATE.out.sample_outputs.map { r -> tuple(r.meta, r.aln) }
    } else if (use_roary) {
        ROARY(gff)
        ch_sample_outputs = ROARY.out.sample_outputs
        ch_aln = ROARY.out.sample_outputs.map { r -> tuple(r.meta, r.aln) }
    } else {
        PANAROO(gff)
        ch_sample_outputs = PANAROO.out.sample_outputs
        ch_aln = PANAROO.out.sample_outputs.map { r -> tuple(r.meta, r.filtered_aln) }
    }

    // Per-sample SNP distances
    ch_unmasked_aln = ch_aln.map({ _meta, aln ->
        tuple([name: "core-genome.distance", process_name: "snpdists"], aln)
    })
    SNPDISTS(ch_unmasked_aln)

    emit:
    sample_outputs = ch_sample_outputs
    snpdists_outputs = SNPDISTS.out.sample_outputs
}
