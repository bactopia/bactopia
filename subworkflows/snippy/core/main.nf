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
 * @output aln              Core SNP alignment in FASTA format (polymorphic sites only)
 * @output full_aln         Full core alignment including monomorphic sites
 * @output clean_full_aln   Cleaned full alignment with constant sites for phylogenetic inference
 * @output tab              Core SNPs in TAB format
 * @output vcf              Core SNPs in VCF format
 * @output txt              Core summary statistics (number of SNPs, core genome size)
 * @output samples          List of samples included in the core alignment
 * @output tsv              Pairwise SNP distance matrix from snp-dists
 * @output results          Aggregated results channel containing all output files
 * @output logs             Aggregated logs channel containing all execution logs
 * @output nf_logs          Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions         Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SNIPPY_CORE as SNIPPY_CORE_MODULE  } from '../../../modules/snippy/core/main'
include { SNPDISTS                           } from '../../snpdists/main'
include { flattenPaths                       } from 'plugin/nf-bactopia'
include { gather                             } from 'plugin/nf-bactopia'

workflow SNIPPY_CORE {
    take:
    alignments: Channel<Map, Set<Path>, Set<Path>>
    reference: Path
    mask: Path?

    main:
    SNIPPY_CORE_MODULE(alignments, reference, mask)

    // Per-sample SNP distances
    ch_unmasked_aln = SNIPPY_CORE_MODULE.out.clean_full_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-snp.distance", process_name: "snpdists"], aln]}
    SNPDISTS(ch_unmasked_aln)

    emit:
    aln = SNIPPY_CORE_MODULE.out.aln
    full_aln = SNIPPY_CORE_MODULE.out.full_aln
    clean_full_aln = SNIPPY_CORE_MODULE.out.clean_full_aln
    tab = SNIPPY_CORE_MODULE.out.tab
    vcf = SNIPPY_CORE_MODULE.out.vcf
    txt = SNIPPY_CORE_MODULE.out.txt
    samples = SNIPPY_CORE_MODULE.out.samples
    tsv = SNPDISTS.out.tsv

    // Generic aggregate outputs
    results = SNIPPY_CORE_MODULE.out.supplemental.mix(
        SNIPPY_CORE_MODULE.out.aln,
        SNIPPY_CORE_MODULE.out.full_aln,
        SNIPPY_CORE_MODULE.out.clean_full_aln,
        SNIPPY_CORE_MODULE.out.tab,
        SNIPPY_CORE_MODULE.out.vcf,
        SNIPPY_CORE_MODULE.out.txt,
        SNIPPY_CORE_MODULE.out.samples,
        SNPDISTS.out.results
    )
    logs = SNIPPY_CORE_MODULE.out.logs.mix(
        SNPDISTS.out.logs
    )
    nf_logs = SNIPPY_CORE_MODULE.out.nf_begin.mix(
        SNIPPY_CORE_MODULE.out.nf_err,
        SNIPPY_CORE_MODULE.out.nf_log,
        SNIPPY_CORE_MODULE.out.nf_out,
        SNIPPY_CORE_MODULE.out.nf_run,
        SNIPPY_CORE_MODULE.out.nf_sh,
        SNIPPY_CORE_MODULE.out.nf_trace,
        SNPDISTS.out.nf_logs
    )
    versions = SNIPPY_CORE_MODULE.out.versions.mix(
        SNPDISTS.out.versions
    )
}
