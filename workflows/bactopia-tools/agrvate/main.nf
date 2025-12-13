#!/usr/bin/env nextflow
/**
 * Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
 *
 * This Bactopia Tool uses [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) to rapidly
 * identify the _agr_ locus type and detect agr operon variants in _Staphylococcus aureus_ assemblies.
 * The agr system is a key quorum-sensing regulator of virulence in S. aureus.
 *
 * @status stable
 * @keywords staphylococcus aureus, agr, typing, virulence, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation agrvate
 *
 * @subworkflows bactopiatool_init, agrvate
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                 Agr locus typing results in TSV format
 * @publish supplemental/*        Supplemental output files including detailed agr analysis
 *
 * @section Merged Results
 * @publish agrvate.tsv            Combined agr typing results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { AGRVATE           } from '../../../subworkflows/agrvate/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    // Execute subworkflows
    BACTOPIATOOL_INIT()
    AGRVATE(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(AGRVATE.out.results)
    ch_logs = ch_logs.mix(AGRVATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(AGRVATE.out.nf_logs)
    ch_versions = ch_versions.mix(AGRVATE.out.versions)

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
