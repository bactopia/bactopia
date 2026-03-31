/**
 * Generate core-genome SNP alignment from per-sample Snippy outputs.
 *
 * This subworkflow aggregates individual Snippy variant calls to produce a core-genome
 * alignment using [snippy-core](https://github.com/tseemann/snippy). It identifies core
 * SNPs present across all samples, generates a clean alignment suitable for phylogenetic
 * analysis, and calculates pairwise SNP distances using snp-dists. The output can be
 * used directly with tree-building tools like IQ-TREE, RAxML, or Gubbins.
 *
 * @status stable
 * @keywords variant calling, core genome, snp, alignment, phylogenetics
 * @tags complexity:moderate input-type:multiple output-type:multiple
 * @citation snippy, snpdists
 *
 * @subworkflows snpdists
 * @modules snippy_core
 *
 * @input alignments
 * Channel containing per-sample aligned FASTA files and VCFs from Snippy runs
 *
 * @input reference
 * Reference genome in GenBank or FASTA format used for variant calling
 *
 * @input mask
 * Optional BED file of regions to mask from the core alignment (e.g., recombinant regions, repeat regions)
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `aln`: Core SNP alignment in FASTA format (polymorphic sites only)
 * - `full_aln`: Full core alignment including monomorphic sites
 * - `clean_full_aln`: Cleaned full alignment with constant sites for phylogenetic inference
 * - `tab`: Core SNPs in TAB format
 * - `vcf`: Core SNPs in VCF format
 * - `txt`: Core summary statistics (number of SNPs, core genome size)
 * - `samples`: List of samples included in the core alignment
 * - `supplemental`: Individual sample alignments and intermediate files
 * - `tsv`: Pairwise SNP distance matrix from snp-dists
 */
nextflow.preview.types = true

include { SNIPPY_CORE as SNIPPY_CORE_MODULE } from '../../../modules/snippy/core/main'
include { SNPDISTS                          } from '../../snpdists/main'

workflow SNIPPY_CORE {
    take:
    alignments: Channel<Record>
    reference: Path
    mask: Path?

    main:
    SNIPPY_CORE_MODULE(alignments, reference, mask)

    // Per-sample SNP distances
    SNPDISTS(SNIPPY_CORE_MODULE.out.map { r ->
        record(meta: [name: 'core-snp.distance', process_name: 'snpdists'], aln: r.clean_full_aln)
    })

    emit: // bactopia-lint: ignore S005, S010
    // Downstream inputs
    alignment = SNIPPY_CORE_MODULE.out.map { r ->
        record(meta: [name: "core-snp", process_name: "iqtree"], aln: r.clean_full_aln)
    }
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = SNIPPY_CORE_MODULE.out.mix(SNPDISTS.out.run_outputs)
}
