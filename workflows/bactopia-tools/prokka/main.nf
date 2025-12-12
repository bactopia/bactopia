#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia Tool: Prokka.
 *
 * Whole genome annotation of small genomes (bacterial, archaeal, viral)
 * The `prokka` module uses [Prokka](https://github.com/tseemann/prokka) to rapidly annotate bacterial
 * genomes in a standardized fashion.
 *
 * @status stable
 *
 * @subworkflows bactopiatool_init, prokka
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @section Per-Sample Results
 * @publish *.gff    Genome annotation in GFF3 format
 * @publish *.gbk    Genome annotation in GenBank format
 * @publish *.faa    Protein sequences
 * @publish *.fna    Nucleotide sequences
 * @publish *.ffn    Feature nucleotide sequences
 *
 * @section Merged Results
 *
 * @publish merged-*    Aggregated results from all samples
 *
 * @section Execution Logs
 *
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 *
 * @publish versions.yml Software version information
   */

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    prokka_proteins    : Path?
    prokka_prodigal_tf : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PROKKA            } from '../../../subworkflows/prokka/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    PROKKA(
        BACTOPIATOOL_INIT.out.samples,
        params.prokka_proteins,
        params.prokka_prodigal_tf
    )

    // Collect outputs
    ch_results = ch_results.mix(PROKKA.out.results)
    ch_logs = ch_logs.mix(PROKKA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PROKKA.out.nf_logs)
    ch_versions = ch_versions.mix(PROKKA.out.versions)

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
