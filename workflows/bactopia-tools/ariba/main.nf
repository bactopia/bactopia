#!/usr/bin/env nextflow
/**
 * Gene identification through local assemblies.
 *
 * This Bactopia Tool uses [ARIBA](https://github.com/sanger-pathogens/ariba) to rapidly
 * identify genes in a database by creating local assemblies from short-read data.
 * ARIBA performs reference-based assembly and variant calling for gene detection.
 *
 * @status stable
 * @keywords fastq, assembly, resistance, virulence, gene detection, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,database-dependent,gene-identification
 * @citation ariba
 *
 * @subworkflows bactopiatool_init, ariba
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input ariba_db
 * Name of the ARIBA database to use for gene identification
 *
 * @section Per-Sample Results
 * @publish *-report.tsv           Gene detection report for each sample
 * @publish *-summary.csv          Summary of gene detection results
 * @publish assembled_genes.fa.gz  Assembled genes in compressed FASTA format
 * @publish assembled_seqs.fa.gz   Assembled sequences matching references
 * @publish assemblies.fa.gz       Raw local assemblies
 * @publish debug.report.tsv       Detailed report including synonymous mutations
 * @publish log.clusters.gz        Analysis log file
 * @publish version_info.txt       Version information for ARIBA and dependencies
 *
 * @section Merged Results
 * @publish ariba-report.tsv       Merged gene detection reports from all samples
 * @publish ariba-summary.csv      Merged summary reports from all samples
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml           Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    ariba_db : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ARIBA             } from '../../../subworkflows/ariba/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ARIBA(
        BACTOPIATOOL_INIT.out.samples,
        params.ariba_db
    )

    // Collect outputs
    ch_results = ch_results.mix(ARIBA.out.results)
    ch_logs = ch_logs.mix(ARIBA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ARIBA.out.nf_logs)
    ch_versions = ch_versions.mix(ARIBA.out.versions)

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
