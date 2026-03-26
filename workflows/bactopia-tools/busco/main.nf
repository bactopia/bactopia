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
 * @subworkflows bactopiatool_init, busco
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
    rundir   : String

    // Tool-specific parameters
    busco_lineage : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { BUSCO             } from '../../../subworkflows/busco/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BUSCO(
        BACTOPIATOOL_INIT.out.assembly,
        params.busco_lineage
    )

    ch_sample_nf_logs = BUSCO.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = BUSCO.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = BUSCO.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = BUSCO.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results      >> "${r.meta.output_dir}/"
            r.supplemental >> "${r.meta.output_dir}/"
            r.logs         >> "${r.meta.logs_dir}/"
            r.versions     >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
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
