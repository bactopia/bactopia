#!/usr/bin/env nextflow
/**
 * Quality assessment of assembled contigs using QUAST.
 *
 * This Bactopia Tool uses [QUAST](https://github.com/ablab/quast) to evaluate the quality
 * of assembled contigs. QUAST (Quality Assessment Tool for Genome Assemblies) generates
 * comprehensive reports including numerous plots and tables that help assess assembly
 * quality metrics such as N50, GC content, genome fraction, and misassembly rates.
 * It produces both per-sample assessments and merged summaries for comparative analysis
 * across multiple samples.
 *
 * @status stable
 * @keywords assembly, quality, assessment, metrics, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, quast
 *
 * @subworkflows bactopiatool_init, quast
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                      Summary statistics of QUAST assessment for each sample
 * @publish basic_stats/               Directory containing plots of assembly metrics (GC content, NGx, Nx)
 * @publish icarus.html                Icarus main menu with links to interactive viewers
 * @publish icarus_viewers/            Additional reports and viewers for Icarus
 * @publish predicted_genes/           Directory containing predicted gene information
 * @publish report.*                   Assessment summary in various formats (html, pdf, tex, tsv, txt)
 * @publish transposed_report.*        Transposed version of the assessment summary (tex, tsv, txt)
 *
 * @section Merged Results
 * @publish quast.tsv                  Merged TSV file with QUAST summary statistics from all samples
 *
 * @section Execution Logs
 * @publish logs/**                    Tool execution logs including QUAST logs
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { QUAST             } from '../../../subworkflows/quast/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    QUAST(BACTOPIATOOL_INIT.out.samples_2)

    // Collect outputs
    ch_results = ch_results.mix(QUAST.out.results)
    ch_logs = ch_logs.mix(QUAST.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QUAST.out.nf_logs)
    ch_versions = ch_versions.mix(QUAST.out.versions)

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
