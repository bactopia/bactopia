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

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { QUAST               } from '../../../subworkflows/quast/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_quast = QUAST(ch_bactopiatool.assembly_meta)

    publish:
    // Per-sample
    sample_outputs = ch_quast.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_quast.sample_outputs)
    // Run-level
    run_outputs = ch_quast.run_outputs
    run_nf_logs = collectNextflowLogs(ch_quast.run_outputs)
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
