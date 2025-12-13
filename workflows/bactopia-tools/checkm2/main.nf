#!/usr/bin/env nextflow
/**
 * Machine learning-based assessment of microbial genome assembly quality.
 *
 * This Bactopia Tool uses [CheckM2](https://github.com/chklovski/CheckM2) to assess the quality
 * of microbial genomes recovered from isolates, single cells, and metagenomes using
 * advanced machine learning approaches.
 *
 * @status stable
 * @keywords assembly quality, microbial genomes, machine learning, completeness, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,assembly-qa,machine-learning
 * @citation checkm2
 *
 * @subworkflows bactopiatool_init, checkm2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input checkm2_db
 * Path to CheckM2 database
 *
 * @input download_checkm2
 * Download CheckM2 database if not found locally
 *
 * @section Quality Assessment
 * @publish diamond_output/**    Directory with intermediate results from CheckM2 processing
 * @publish protein_files/**     Directory containing protein files used for analysis
 * @publish quality_report.tsv   Output file with completeness statistics
 *
 * @section Merged Results
 * @publish checkm2.tsv          Merged TSV file with CheckM2 results from all samples
 *
 * @section Execution Logs
 * @publish logs/**              Tool execution logs
 * @publish logs/nf-*            Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml         Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    checkm2_db        : Path
    download_checkm2  : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { CHECKM2           } from '../../../subworkflows/checkm2/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    CHECKM2(
        BACTOPIATOOL_INIT.out.samples,
        params.checkm2_db,
        params.download_checkm2
    )

    // Collect outputs
    ch_results = ch_results.mix(CHECKM2.out.results)
    ch_logs = ch_logs.mix(CHECKM2.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CHECKM2.out.nf_logs)
    ch_versions = ch_versions.mix(CHECKM2.out.versions)

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
