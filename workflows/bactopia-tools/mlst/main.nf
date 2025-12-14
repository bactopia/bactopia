#!/usr/bin/env nextflow
/**
 * Automatic Multi-Locus Sequence Type (MLST) calling from assembled contigs.
 *
 * This Bactopia Tool uses [mlst](https://github.com/tseemann/mlst) to scan genome assemblies
 * and determine the sequence type based on [PubMLST](https://pubmlst.org/) schemes. The workflow
 * automatically detects the appropriate MLST scheme for each organism and provides
 * standardized sequence type assignments.
 *
 * @status stable
 * @keywords typing, MLST, sequence type, alleles, PubMLST
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, mlst, pubmlst
 *
 * @subworkflows bactopiatool_init, mlst
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                      Tab-delimited file with MLST results including scheme, ST, and allele profiles
 *
 * @section Merged Results
 * @publish mlst.tsv                   Merged TSV file containing MLST results from all samples
 *
 * @section Execution Logs
 * @publish logs/mlst/*                Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    mlst_db : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MLST              } from '../../../subworkflows/mlst/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    
    // MLST database can be optionally provided via params
    MLST(
        BACTOPIATOOL_INIT.out.assembly,
        params.mlst_db
    )

    // Collect outputs
    ch_results = ch_results.mix(MLST.out.results)
    ch_logs = ch_logs.mix(MLST.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MLST.out.nf_logs)
    ch_versions = ch_versions.mix(MLST.out.versions)

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
