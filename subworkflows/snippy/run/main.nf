/**
 * Call variants against a reference genome using Snippy.
 *
 * This subworkflow performs rapid haploid variant calling from bacterial sequence reads
 * using [Snippy](https://github.com/tseemann/snippy). It maps reads to a reference genome,
 * identifies SNPs and indels, and generates consensus sequences. The tool produces multiple
 * output formats including VCF, aligned FASTA, and annotated variants for downstream
 * phylogenetic analysis with snippy-core.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords variant calling, snp, reference mapping, phylogenetics, outbreak
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation snippy
 *
 * @modules snippy_run
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input reference
 * Reference genome in GenBank format (preferred, for annotation) or FASTA format
 *
 * @output sample_outputs
 *   - `aligned_fa`: A version of the reference with - at zero coverage positions
 *   - `vcf`: The final annotated variants in VCF format
 *   - `aligned_fa_error`: Aligned FASTA file generated during error state
 *   - `vcf_error`: VCF file generated during error state
 *   - `error`: Error log text file
 *   - `annotated_vcf`: Annotated VCF file
 *   - `bam`: The alignments in BAM format (includes unmapped/multimapping)
 *   - `bai`: Index for the BAM file
 *   - `bed`: The variants in BED format
 *   - `consensus_fa`: Reference genome with all variants instantiated
 *   - `consensus_subs_fa`: Reference genome with only substitution variants instantiated
 *   - `consensus_subs_masked_fa`: Reference genome with substitutions instantiated and low coverage masked
 *   - `coverage`: Per-base coverage depth information
 *   - `csv`: A comma-separated summary of variants
 *   - `filt_vcf`: The filtered variant calls from Freebayes
 *   - `gff`: The variants in GFF3 format
 *   - `html`: A HTML summary of the variants
 *   - `raw_vcf`: The unfiltered variant calls from Freebayes
 *   - `subs_vcf`: VCF containing only substitution variants
 *   - `tab`: A simple tab-separated summary of all variants
 *   - `txt`: Tab-separated columnar list of alignment statistics
 */
nextflow.preview.types = true

include { SNIPPY_RUN } from '../../../modules/snippy/run/main'
include { gather     } from 'plugin/nf-bactopia'

workflow SNIPPY {
    take:
    reads : Channel<Record>
    reference : Path

    main:
    SNIPPY_RUN(reads, reference)

    emit:
    sample_outputs = SNIPPY_RUN.out
}
