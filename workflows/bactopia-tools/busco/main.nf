#!/usr/bin/env nextflow
/**
 * Assessment of genome assembly completeness using evolutionarily informed expectations.
 *
 * This Bactopia Tool uses [BUSCO](https://gitlab.com/ezlab/busco) (Benchmarking Universal Single-Copy Orthologs)
 * to assess the completeness of genome assemblies by searching for single-copy orthologs. The workflow
 * processes each assembly against a specified lineage dataset and provides comprehensive completeness metrics.
 *
 * @status stable
 * @keywords assembly, completeness, assessment, orthologs, quality control
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation,database-dependent
 * @citation busco, csvtk
 *
 * @subworkflows bactopiatool_init, busco
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input busco_lineage
 * BUSCO lineage dataset for completeness assessment
 *
 * @section Per-Sample Results
 * @publish run_                      BUSCO analysis output directory for each lineage
 * @publish run_/full_table.tsv        Complete results with scores and lengths of BUSCO matches
 * @publish run_/missing_busco_list.tsv List of missing BUSCO genes
 * @publish run_/short_summary.txt     Summary of BUSCO assessment results
 * @publish run_/short_summary.json    Summary of BUSCO assessment in JSON format
 * @publish *-summary.txt             Per-sample BUSCO summary file
 * @publish *-summary.json            Per-sample BUSCO summary in JSON format
 *
 * @section Merged Results
 * @publish busco.tsv                   Merged TSV file containing BUSCO summaries from all samples
 *
 * @section Execution Logs
 * @publish logs/busco/*                Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                   Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    busco_lineage : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { BUSCO             } from '../../../subworkflows/busco/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    BUSCO(
        BACTOPIATOOL_INIT.out.samples,
        params.busco_lineage
    )

    // Collect outputs
    ch_results = ch_results.mix(BUSCO.out.results)
    ch_logs = ch_logs.mix(BUSCO.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BUSCO.out.nf_logs)
    ch_versions = ch_versions.mix(BUSCO.out.versions)

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
