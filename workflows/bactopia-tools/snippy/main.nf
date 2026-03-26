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
 * @subworkflows bactopiatool_init, ncbigenomedownload, snippy, snippy_core, gubbins, iqtree
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
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    reference          : Path
    accession          : String
    snippy_core_mask   : Path
    skip_recombination : Boolean
    skip_phylogeny     : Boolean
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { SNIPPY             } from '../../../subworkflows/snippy/run/main'
include { SNIPPY_CORE        } from '../../../subworkflows/snippy/core/main'
include { GUBBINS            } from '../../../subworkflows/gubbins/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'

workflow {
    main:
    BACTOPIATOOL_INIT()

    // Download if applicable
    ch_reference = channel.empty()
    if (params.reference) {
        ch_reference = [[id: 'snippy'], params.reference]
    } else if (params.accession) {
        NCBIGENOMEDOWNLOAD([])
        ch_reference = NCBIGENOMEDOWNLOAD.out.bactopia_tools.first()
    }

    // Run Snippy per-sample
    SNIPPY(BACTOPIATOOL_INIT.out.reads, ch_reference)
    ch_sample_outputs = SNIPPY.out.sample_outputs

    // Collect VCFs and aligned FAs across all samples for core-SNP analysis
    ch_vcfs = SNIPPY.out.sample_outputs.flatMap { r -> r.vcf }.collect()
    ch_fas = SNIPPY.out.sample_outputs.flatMap { r -> r.aligned_fa }.collect()
    ch_snippy_core_input = ch_vcfs.combine(ch_fas).map { vcfs, fas ->
        record(_meta: [id: 'core-snp'], _vcf: vcfs.toSet(), _aligned_fa: fas.toSet())
    }

    // Identify core SNPs
    SNIPPY_CORE(ch_snippy_core_input, ch_reference, params.snippy_core_mask ? [params.snippy_core_mask] : [])
    ch_sample_outputs = ch_sample_outputs
        .mix(SNIPPY_CORE.out.sample_outputs)
        .mix(SNIPPY_CORE.out.snpdists_outputs)

    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        ch_gubbins_input = SNIPPY_CORE.out.sample_outputs.map { r ->
            record(_meta: r.meta, msa: r.clean_full_aln)
        }
        GUBBINS(ch_gubbins_input)
        ch_sample_outputs = ch_sample_outputs
            .mix(GUBBINS.out.sample_outputs)
            .mix(GUBBINS.out.snpdists_outputs)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        ch_final_aln = channel.empty()
        if (!params.skip_recombination) {
            ch_final_aln = GUBBINS.out.sample_outputs.map { r ->
                record(_meta: [name: "core-snp", process_name: "iqtree"], _msa: [r.masked_aln].toSet())
            }
        } else {
            ch_final_aln = SNIPPY_CORE.out.sample_outputs.map { r ->
                record(_meta: [name: "core-snp", process_name: "iqtree"], _msa: [r.clean_full_aln].toSet())
            }
        }
        IQTREE(ch_final_aln)
        ch_sample_outputs = ch_sample_outputs.mix(IQTREE.out.sample_outputs)
    }

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ch_sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ch_sample_outputs
    sample_nf_logs = ch_sample_nf_logs
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
}
