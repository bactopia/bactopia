#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia Tool: Pangenome.
 *
 * Pangenome analysis with optional core-genome phylogeny
 * The `pangenome` subworkflow allows you to create a pan-genome with [PIRATE](https://github.com/SionBayliss/PIRATE),
 * [Panaroo](https://github.com/gtonkinhill/panaroo), or [Roary](https://github.com/sanger-pathogens/Roary)) of your samples.
 * You can further supplement your pan-genome by including completed genomes. This is possible using the `--species`
 * or `--accessions` parameters. If used, [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) will
 * download available completed genomes available from RefSeq. Any downloaded genomes will be annotated with
 * [Prokka](https://github.com/tseemann/prokka) to create compatible GFF3 files.
 * A phylogeny, based on the core-genome alignment, will be created by [IQ-Tree](https://github.com/Cibiv/IQ-TREE). Optionally
 * a recombination-masked core-genome alignment can be created with [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML)
 * and [maskrc-svg](https://github.com/kwongj/maskrc-svg).
 * Finally, the core genome pair-wise SNP distance for each sample is also calculated with
 * [snp-dists](https://github.com/tseemann/snp-dists) and additional pan-genome wide association studies can be conducted
 * using [Scoary](https://github.com/AdmiralenOla/Scoary).
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny
 *
 * @subworkflows bactopiatool_init, pangenome
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @section Per-Sample Results
 * @publish *    Analysis results
 *
 * @section Merged Results
 * @publish merged-*    Aggregated results from all samples
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution logs
 *
 * @section Versions
 * @publish versions.yml Software version information
   */

params {
    rundir : String
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { PROKKA             } from '../../../subworkflows/prokka/main'
include { PANGENOME          } from '../../../subworkflows/pangenome/main'
include { CLONALFRAMEML      } from '../../../subworkflows/clonalframeml/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'
include { SCOARY             } from '../../../subworkflows/scoary/main'
include { formatSamples      } from 'plugin/nf-bactopia'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ch_samples = BACTOPIATOOL_INIT.out.samples

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
