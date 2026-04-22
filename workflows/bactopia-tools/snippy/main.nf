#!/usr/bin/env nextflow
/**
 * Rapid haplotype variant calling and core genome alignment.
 *
 * This Bactopia Tool uses [Snippy](https://github.com/tseemann/snippy) to find SNPs between a
 * reference genome and a set of reads, perform core genome alignment, and generate
 * phylogenetic trees. It includes optional recombination detection with Gubbins
 * and phylogenetic tree construction with IQ-Tree.
 *
 * @status stable
 * @keywords snp, variant calling, phylogeny, core genome, snippy, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,comparative,phylogeny
 * @citation snippy, gubbins, iqtree
 *
 * @subworkflows utils_bactopia-tools, ncbigenomedownload, snippy_run, snippy_core, gubbins, iqtree
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input reference
 * Path to reference FASTA file for variant calling
 *
 * @input accession
 * NCBI Assembly RefSeq accession to use as reference
 *
 * @input snippy_core_mask
 * Path to BED file containing core genome regions
 *
 * @input skip_recombination
 * Skip recombination detection with Gubbins
 *
 * @input skip_phylogeny
 * Skip phylogenetic tree construction
 *
 * @section Variant Calling
 * @publish *.vcf            Variant calls in VCF format
 * @publish *.bam            Alignment file
 * @publish *.txt            Snippy summary report
 *
 * @section Core Genome Alignment
 * @publish core.full.aln     Full core genome alignment
 * @publish core.snps.aln     Core SNP alignment
 *
 * @section Recombination Analysis
 * @note Only created if recombination analysis is enabled
 * @publish *.filtered.aln    Alignment with recombination regions removed
 * @publish *.gff            Recombination predictions
 *
 * @section Phylogeny
 * @note Only created if phylogeny analysis is enabled
 * @publish *.treefile       Phylogenetic tree in Newick format
 *
 * @section Merged Results
 * @publish snippy.tsv        Merged summary of Snippy analyses
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    reference          : Path?
    accession          : String
    snippy_core_mask   : Path?
    skip_recombination : Boolean
    skip_phylogeny     : Boolean
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { NCBIGENOMEDOWNLOAD  } from '../../../subworkflows/ncbigenomedownload/main'
include { SNIPPY              } from '../../../subworkflows/snippy/run/main'
include { SNIPPY_CORE         } from '../../../subworkflows/snippy/core/main'
include { GUBBINS             } from '../../../subworkflows/gubbins/main'
include { IQTREE              } from '../../../subworkflows/iqtree/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'
include { gatherFields        } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()

    // Download if applicable
    ch_reference = null
    if (params.reference) {
        ch_reference = params.reference
    } else if (params.accession) {
        ch_ncbigenomedownload = NCBIGENOMEDOWNLOAD(null)
        ch_reference = ch_ncbigenomedownload.reference
    }

    // Run Snippy per-sample
    ch_snippy = SNIPPY(ch_bactopiatool.reads, ch_reference)
    ch_sample_outputs = ch_snippy.sample_outputs
    ch_run_outputs = ch_snippy.run_outputs

    // Collect per-sample VCFs and aligned FAs for core-SNP analysis
    ch_core_input = gatherFields(
        ch_snippy.variants,
        [vcf: '_vcf', aligned_fa: '_aligned_fa'],
        [name: 'core-snp']
    )

    // Identify core SNPs
    ch_snippy_core = SNIPPY_CORE(ch_core_input, ch_reference, params.snippy_core_mask)
    ch_sample_outputs = ch_sample_outputs.mix(ch_snippy_core.sample_outputs)
    ch_run_outputs = ch_run_outputs.mix(ch_snippy_core.run_outputs)

    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        ch_gubbins = GUBBINS(ch_snippy_core.alignment)
        ch_sample_outputs = ch_sample_outputs.mix(ch_gubbins.sample_outputs)
        ch_run_outputs = ch_run_outputs.mix(ch_gubbins.run_outputs)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        ch_iqtree = IQTREE(params.skip_recombination ? ch_snippy_core.alignment : ch_gubbins.alignment)
        ch_sample_outputs = ch_sample_outputs.mix(ch_iqtree.sample_outputs)
        ch_run_outputs = ch_run_outputs.mix(ch_iqtree.run_outputs)
    }

    publish:
    // Per-sample
    sample_outputs = ch_sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_sample_outputs)
    // Run-level
    run_outputs = ch_run_outputs
    run_nf_logs = collectNextflowLogs(ch_run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
