#!/usr/bin/env nextflow
/**
 * Taxonomic classification of metagenomic sequence reads.
 *
 * This Bactopia Tool uses [Kraken2](https://github.com/DerrickWood/kraken2) to assign taxonomic
 * classifications to metagenomic sequence reads. It creates reports of taxonomic assignments
 * and generates interactive Krona plots for visualization of classification results.
 *
 * @status stable
 * @keywords taxonomy, metagenomics, classification, fastq, kraken2, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,classification,visualization
 * @citation kraken2
 *
 * @subworkflows bactopiatool_init, kraken2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input kraken2_db
 * Path to Kraken2 database for taxonomic classification
 *
 * @section Per-Sample Results
 * @publish *.report          Kraken2 taxonomic classification report
 * @publish *.kraken2         Kraken2 classification results
 * @publish *.html            Krona interactive HTML report
 *
 * @section Merged Results
 * @publish merged-report     Merged Kraken2 reports from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    kraken2_db : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { KRAKEN2           } from '../../../subworkflows/kraken2/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    KRAKEN2(
        BACTOPIATOOL_INIT.out.samples,
        params.kraken2_db
    )

    // Collect outputs
    ch_results = ch_results.mix(KRAKEN2.out.results)
    ch_logs = ch_logs.mix(KRAKEN2.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KRAKEN2.out.nf_logs)
    ch_versions = ch_versions.mix(KRAKEN2.out.versions)

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
