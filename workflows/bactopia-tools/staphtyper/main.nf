#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia Tool: Staphtyper.
 *
 * Determine the agr, spa and SCCmec types for _Staphylococcus aureus_ genomes
 * The `staphtyper` subworkflow includes multiple tools that are specific for typing certain features
 * of *Staphylococcus aureus*. Currently `staphtyper` includes
 * 1. [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) - *agr* locus type and *agr* operon variants.
 * 2. [spaTyper](https://github.com/HCGB-IGTP/spaTyper) - *spa* type
 * 3. [sccmec](https://github.com/rpetit3/sccmec) - SCCmec type
 * This tool will evolve with *S. aureus* genomics, so you can expect it to add more typing methods
 * (maybe even replace current methods) in the future. If a certain typing method for *S. aureus*
 * please feel free to suggest it be added!~
 *
 * @status stable
 *
 * @subworkflows bactopiatool_init, staphtyper
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
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml Software version information
   */

params {
    rundir   : String

    // Tool-specific parameters
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { STAPHTYPER        } from '../../../subworkflows/staphtyper/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    STAPHTYPER(
        BACTOPIATOOL_INIT.out.samples,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )

    // Collect outputs
    ch_results = ch_results.mix(STAPHTYPER.out.results)
    ch_logs = ch_logs.mix(STAPHTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(STAPHTYPER.out.nf_logs)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)

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
