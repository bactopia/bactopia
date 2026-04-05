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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,conditional-logic
 * @citation pirate, panaroo, roary, snpdists
 *
 * @subworkflows pirate, roary, panaroo, snpdists
 *
 * @input record(meta, gff)
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
 *
 * @output run_outputs
 * - `aln`: Core-genome alignment in FASTA format
 * - `csv`: Gene presence/absence matrix
 * - `supplemental`: Intermediate files and detailed outputs
 * - `tsv`: Pairwise SNP distance matrix from core-genome alignment
 *
 * @output alignment
 * - `aln`: Core-genome alignment for downstream analysis (e.g., recombination detection)
 *
 * @output phylogeny_input
 * - `aln`: Core-genome alignment with iqtree-ready meta for phylogeny construction
 *
 * @output csv
 * - `csv`: Gene presence/absence matrix for downstream analysis (e.g., pan-GWAS)
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
    ch_run_outputs = channel.empty()

    // Choose pangenome tool based on params
    if (use_pirate) {
        PIRATE(gff)
        ch_run_outputs = PIRATE.out.run_outputs
    } else if (use_roary) {
        ROARY(gff)
        ch_run_outputs = ROARY.out.run_outputs
    } else {
        PANAROO(gff)
        ch_run_outputs = PANAROO.out.run_outputs
    }

    // SNP distances (panaroo uses filtered_aln, others use aln)
    SNPDISTS(ch_run_outputs.map { r ->
        def core_aln = use_pirate || use_roary ? r.aln : r.filtered_aln
        record(meta: [name: 'core-genome.distance', process_name: 'snpdists'], aln: core_aln)
    })

    emit: // bactopia-lint: ignore S005, S010
    // Downstream inputs
    alignment = ch_run_outputs.map { r ->
        record(meta: r.meta, aln: (use_pirate || use_roary ? r.aln : r.filtered_aln))
    }
    phylogeny_input = ch_run_outputs.map { r ->
        record(meta: [name: "core-genome", process_name: "iqtree"], aln: (use_pirate || use_roary ? r.aln : r.filtered_aln))
    }
    csv = ch_run_outputs.map { r -> record(meta: r.meta, csv: r.csv) }
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_run_outputs.mix(SNPDISTS.out.run_outputs)
}
