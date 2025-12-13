#!/usr/bin/env nextflow
/**
 * Mass screening of contigs for antimicrobial resistance and virulence genes.
 *
 * This Bactopia Tool uses [Abricate](https://github.com/tseemann/abricate) to screen
 * assemblies against multiple resistance and virulence gene databases, including
 * NCBI, CARD, RESFINDER, ARG-ANNOT, VFDB, PLASMIDFINDER, ECOLI_VF, and MEGARES.
 * It processes a Bactopia analysis directory, runs Abricate on each sample, and
 * creates a merged summary report.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, virulence, screening, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation abricate, arg_annot, card, csvtk, ecoh, megares2, ncbi_reference_gene_catalog, plasmidfinder, resfinder, vfdb
 *
 * @subworkflows bactopiatool_init, abricate
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv            Tab-delimited report of Abricate screening results for each sample
 *
 * @section Merged Results
 * @publish abricate.tsv     Merged TSV report containing Abricate results from all samples
 *
 * @section Execution Logs
 * @publish logs/abricate/*  Tool execution logs (stdout/stderr)
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ABRICATE          } from '../../../subworkflows/abricate/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ABRICATE(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(ABRICATE.out.results)
    ch_logs = ch_logs.mix(ABRICATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ABRICATE.out.nf_logs)
    ch_versions = ch_versions.mix(ABRICATE.out.versions)

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
