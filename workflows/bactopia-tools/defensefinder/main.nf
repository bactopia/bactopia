#!/usr/bin/env nextflow
/**
 * Systematic identification of anti-phage defense systems.
 *
 * This Bactopia Tool uses [DefenseFinder](https://github.com/mdmparis/defense-finder)
 * to systematically search for and identify all known anti-phage defense systems
 * in bacterial genomes using HMM-based protein domain detection.
 *
 * @status stable
 * @keywords anti-phage, defense systems, hmm, protein domains, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,hmm-search
 *
 * @subworkflows bactopiatool_init, defensefinder
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.prt                          FASTA file containing all proteins found in defense systems
 * @publish *.prt.idx                      Index file for the proteins file
 * @publish *defense_finder_genes.tsv      TSV file with each gene found in defense systems
 * @publish *defense_finder_hmmer.tsv      TSV file with each HMM hit
 * @publish *defense_finder_systems.tsv    TSV file with information about each system found
 * @publish *.macsydata.tar.gz             Raw MACSyFinder output file (requires --df_preserveraw)
 *
 * @section Merged Results
 * @publish defensefinder-genes.tsv        Merged TSV of all genes found in defense systems
 * @publish defensefinder-hmmer.tsv        Merged TSV of all HMM hits
 * @publish defensefinder-systems.tsv      Merged TSV of all information about systems found
 *
 * @section Execution Logs
 * @publish logs/**                        Tool execution logs
 * @publish logs/nf-*                      Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                   Software version information
 *
 * @citation defensefinder
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { DEFENSEFINDER     } from '../../../subworkflows/defensefinder/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    DEFENSEFINDER(BACTOPIATOOL_INIT.out.assembly)

    // Collect outputs
    ch_results = ch_results.mix(DEFENSEFINDER.out.results)
    ch_logs = ch_logs.mix(DEFENSEFINDER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(DEFENSEFINDER.out.nf_logs)
    ch_versions = ch_versions.mix(DEFENSEFINDER.out.versions)

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
