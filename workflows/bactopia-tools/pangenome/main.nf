#!/usr/bin/env nextflow
/**
 * Pangenome analysis with optional core-genome phylogeny.
 *
 * This Bactopia Tool creates a pangenome from GFF3 annotation files using one of three
 * tools: [Panaroo](https://github.com/gtonkinhill/panaroo) (default),
 * [PIRATE](https://github.com/SionBayliss/PIRATE), or
 * [Roary](https://github.com/sanger-pathogens/roary). It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations.
 * You can supplement your pangenome with completed genomes using the --species or
 * --accessions parameters, which downloads genomes from RefSeq and annotates them with
 * Prokka. A phylogeny based on the core-genome alignment is created by IQ-Tree, with
 * optional recombination masking using ClonalFrameML. Finally, pan-genome wide
 * association studies can be conducted using Scoary.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,aggregation,conditional-logic
 * @citation clonalframeml, iqtree, iqtree_modelfinder, iqtree_ufboot, ncbigenomedownload, panaroo, pirate, prokka, roary, scoary
 *
 * @subworkflows bactopiatool_init, pangenome, ncbigenomedownload, prokka, clonalframeml, iqtree, scoary
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Pangenome Results
 * @publish *.aln                 Core-genome alignment file containing genes present across all input genomes
 * @publish *.csv                 Gene presence/absence matrix showing which genes are present in each genome
 * @publish *.tsv                 SNP distance matrix between all samples
 *
 * @section Phylogeny Results
 * @note Only created if --skip_phylogeny is not enabled
 * @publish *.treefile            Maximum likelihood phylogenetic tree in Newick format
 * @publish *.iqtree              IQ-Tree analysis report with model selection and support values
 * @publish *.log                 IQ-Tree execution log
 *
 * @section Recombination Analysis
 * @note Only created if --skip_recombination is not enabled
 * @publish *.masked.aln          Core-genome alignment with recombination regions masked
 *
 * @section Association Analysis
 * @note Only created if --scoary_traits is specified
 * @publish scoary/*              Scoary association analysis results and plots
 *
 * @section Panaroo Results
 * @note Only created when Panaroo is selected as the pangenome tool
 * @publish panaroo/*             Panaroo-specific output files including graph and statistics
 *
 * @section PIRATE Results
 * @note Only created when PIRATE is selected as the pangenome tool
 * @publish pirate/*              PIRATE-specific output files including gene families and clusters
 *
 * @section Roary Results
 * @note Only created when Roary is selected as the pangenome tool
 * @publish roary/*               Roary-specific output files including gene presence/absence matrices
 *
 * @section Execution Logs
 * @publish logs/pangenome/*      Pangenome tool execution logs (stdout/stderr)
 * @publish logs/clonalframeml/*  ClonalFrameML execution logs (if executed)
 * @publish logs/iqtree/*         IQ-Tree execution logs (if executed)
 * @publish logs/scoary/*         Scoary execution logs (if executed)
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Pangenome tool selection
    use_pirate : Boolean
    use_roary  : Boolean

    // Reference genome parameters
    species   : String?
    accession : String?
    accessions : Path?

    // Prokka parameters
    prokka_proteins    : Path?
    prokka_prodigal_tf : Path?

    // Analysis options
    skip_recombination : Boolean
    skip_phylogeny     : Boolean
    scoary_traits      : Path?
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { PROKKA             } from '../../../subworkflows/prokka/main'
include { PANGENOME          } from '../../../subworkflows/pangenome/main'
include { CLONALFRAMEML      } from '../../../subworkflows/clonalframeml/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'
include { SCOARY             } from '../../../subworkflows/scoary/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ch_samples = BACTOPIATOOL_INIT.out.gffs

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        PROKKA(
            NCBIGENOMEDOWNLOAD.out.bactopia_tools,
            params.prokka_proteins ? file(params.prokka_proteins, checkIfExists: true) : [],
            params.prokka_prodigal_tf ? file(params.prokka_prodigal_tf, checkIfExists: true) : []
        )
        ch_samples = ch_samples.mix(PROKKA.out.gff)
    }

    // Create the pangenome
    ch_merge_gff = ch_samples.collect{_meta, gff -> gff}.map{ gff -> [[id: params.use_pirate? 'pirate' : (params.use_roary ? 'roary' : 'panaroo')], gff]}
    PANGENOME(ch_merge_gff, params.use_pirate, params.use_roary)
    ch_results = ch_results.mix(PANGENOME.out.results)
    ch_logs = ch_logs.mix(PANGENOME.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PANGENOME.out.nf_logs)
    ch_versions = ch_versions.mix(PANGENOME.out.versions)
    
    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        // Run ClonalFrameML
        CLONALFRAMEML(PANGENOME.out.aln)
        ch_results = ch_results.mix(CLONALFRAMEML.out.results)
        ch_logs = ch_logs.mix(CLONALFRAMEML.out.logs)
        ch_nf_logs = ch_nf_logs.mix(CLONALFRAMEML.out.nf_logs)
        ch_versions = ch_versions.mix(CLONALFRAMEML.out.versions)
    }

    // (optional) Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        ch_final_aln = channel.empty()
        if (params.skip_recombination) {
            ch_final_aln = PANGENOME.out.aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome", process_name: "iqtree"], aln]}
        } else {
            ch_final_aln = CLONALFRAMEML.out.masked_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome", process_name: "iqtree"], aln]}
        }
        IQTREE(ch_final_aln)
        ch_results = ch_results.mix(IQTREE.out.results)
        ch_logs = ch_logs.mix(IQTREE.out.logs)
        ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    // Pan-genome GWAS
    if (params.scoary_traits) {
        SCOARY(PANGENOME.out.csv, file(params.scoary_traits, checkIfExists: true))
        ch_results = ch_results.mix(SCOARY.out.csv)
        ch_logs = ch_logs.mix(SCOARY.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SCOARY.out.nf_logs)
        ch_versions = ch_versions.mix(SCOARY.out.versions)
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
