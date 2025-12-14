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
include { formatSamples      } from 'plugin/nf-bactopia'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
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
    SNIPPY(
        BACTOPIATOOL_INIT.out.reads,
        ch_reference
    )
    ch_results = ch_results.mix(SNIPPY.out.results)
    ch_logs = ch_logs.mix(SNIPPY.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNIPPY.out.nf_logs)
    ch_versions = ch_versions.mix(SNIPPY.out.versions)

    // Identify core SNPs
    ch_merge_vcf = SNIPPY.out.vcf.collect{_meta, vcf -> vcf}.map{ vcf -> [[id:'core-snp'], vcf]}
    ch_merge_aligned_fa = SNIPPY.out.aligned_fa.collect{_meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'core-snp'], aligned_fa]}
    ch_snippy_core = ch_merge_vcf.join( ch_merge_aligned_fa )

    // Identify core SNPs
    SNIPPY_CORE(ch_snippy_core, ch_reference, params.snippy_core_mask ? [params.snippy_core_mask] : [])
    ch_results = ch_results.mix(SNIPPY_CORE.out.results)
    ch_logs = ch_logs.mix(SNIPPY_CORE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNIPPY_CORE.out.nf_logs)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)

    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(SNIPPY_CORE.out.clean_full_aln)
        ch_results = ch_results.mix(GUBBINS.out.results)
        ch_logs = ch_logs.mix(GUBBINS.out.logs)
        ch_nf_logs = ch_nf_logs.mix(GUBBINS.out.nf_logs)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE.out.clean_full_aln)
        }
        ch_results = ch_results.mix(IQTREE.out.results)
        ch_logs = ch_logs.mix(IQTREE.out.logs)
        ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
