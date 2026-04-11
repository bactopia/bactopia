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
 * @subworkflows utils_bactopia-tools, busco
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
    rundir : String

    // Tool-specific parameters
    busco_lineage : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BUSCO               } from '../../../subworkflows/busco/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_busco = BUSCO(ch_bactopiatool.assembly, params.busco_lineage)

    publish:
    // Per-sample
    sample_outputs = ch_busco.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_busco.sample_outputs)
    // Run-level
    run_outputs = ch_busco.run_outputs
    run_nf_logs = collectNextflowLogs(ch_busco.run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
